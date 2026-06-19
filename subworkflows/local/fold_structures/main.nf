//
// FOLD_STRUCTURES — group rfnorm XML files by cell_line, optionally run rf-jackknife for
// slope/intercept calibration, run rf-fold to predict secondary structures, then convert
// dot-plot outputs to .bp arc files (genome- and transcript-coordinate) and merge them.
//

include { RNAFRAMEWORK_RFJACKKNIFE                                     } from '../../../modules/local/rnaframework/jackknife/main'
include { RNAFRAMEWORK_RFFOLD                                          } from '../../../modules/local/rnaframework/fold/main'
include { RNAFRAMEWORK_DOTPLOT2BP                                      } from '../../../modules/local/dotplot2bp/main'
include { RNAFRAMEWORK_DOTPLOT2BP as RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT } from '../../../modules/local/dotplot2bp/main'
include { MERGE_BP                                                     } from '../../../modules/local/merge_bp/main'
include { MERGE_BP as MERGE_BP_TRANSCRIPT                              } from '../../../modules/local/merge_bp/main'
include { resolveReferenceKey } from '../../../workflows/rnastructurome_functions.nf'

workflow FOLD_STRUCTURES {

    take:
    ch_rfnorm_xml        // channel: [ val(meta), path(xml) ]  — one entry per rfnorm group
    ch_reference_gtf_map // value:   map of reference_key -> [ val(meta), path(gtf) ]
    pipeline_config      // map

    main:
    ch_versions = channel.empty()

    //
    // Group rfnorm XMLs by cell_line so biological replicates are folded together.
    // The group key intentionally excludes replicate.
    //
    def ch_fold_input = ch_rfnorm_xml
        .map { meta, xml ->
            if (!meta.cell_line) {
                error("Missing cell_line for sample '${meta.id}'. rf-fold replicate grouping requires cell_line.")
            }
            def fold_group = meta.cell_line.toString()
            [ fold_group, [ meta, xml ] ]
        }
        .groupTuple()
        .map { fold_group, entries ->
            def metas      = entries.collect { entry -> entry[0] }
            def xmls       = entries.collect { entry -> entry[1] }.flatten()
            def base       = metas[0]
            def replicates = metas.collect { meta -> (meta.replicate ?: 'na').toString() }.unique().sort()
            def sampleIds  = metas.collect { meta -> (meta.id ?: 'na').toString() }.unique().sort()
            def foldMeta   = base + [
                id              : fold_group,
                fold_group      : fold_group,
                fold_replicates : replicates.join(','),
                fold_source_ids : sampleIds.join(','),
                fold_xml_count  : xmls.size()
            ]
            [ foldMeta, xmls ]
        }

    //
    // Optional rf-jackknife — quality assessment against a reference structure set.
    // When --jackknife_reference is provided, jackknife runs and rf-fold is gated on completion.
    // When --rfjackknife_pool_all is true, the optimal slope/intercept from the jackknife CSV
    // is parsed and injected directly into rf-fold meta, bypassing a separate calibration run.
    //
    def ch_fold_for_rffold = ch_fold_input
    def ch_jackknife_csv   = channel.empty()

    if (pipeline_config.jackknife_reference) {
        def ch_jackknife_reference = channel.value(file(pipeline_config.jackknife_reference.toString(), checkIfExists: true))

        // Jackknife runs per rfnorm group (cell_line + replicate), not per fold group.
        // Fold groups flatten XMLs from multiple replicates which produce identically-named
        // files (e.g. 16S_rRNA.xml); using input*/* staging gives each file its own
        // numbered directory so rf-jackknife receives them as separate experiment dirs.
        def ch_jackknife_input
        if (pipeline_config.rfjackknife_pool_all as Boolean) {
            ch_jackknife_input = ch_rfnorm_xml
                .flatMap { _meta, xmls -> [xmls].flatten() }
                .collect()
                .map { all_xmls -> [ [ id: 'all_groups', fold_group: 'all_groups' ], all_xmls ] }
        } else {
            ch_jackknife_input = ch_rfnorm_xml
        }

        RNAFRAMEWORK_RFJACKKNIFE (
            ch_jackknife_input,
            ch_jackknife_reference
        )
        ch_versions     = ch_versions.mix(RNAFRAMEWORK_RFJACKKNIFE.out.versions.first())
        ch_jackknife_csv = RNAFRAMEWORK_RFJACKKNIFE.out.csv

        if (pipeline_config.rfjackknife_pool_all as Boolean) {
            // Parse optimal slope/intercept from the pooled jackknife CSV and inject into fold meta.
            // FMI.csv is a semicolon-delimited matrix: rows = slopes, columns = intercepts, cells = FMI.
            // Header row: FMI;<intercept0>;<intercept1>;...
            // Data rows:  <slope>;<fmi0>;<fmi1>;...
            def ch_calibration = RNAFRAMEWORK_RFJACKKNIFE.out.csv
                .map { _meta, csv_files ->
                    def csv_file = [csv_files].flatten()[0]
                    def lines = csv_file.readLines()
                    def intercepts = lines[0].split(';').drop(1)*.trim()
                    def best_slope = null; def best_intercept = null; def best_fmi = -1d
                    lines.drop(1).findAll { line -> line.trim() }.each { line ->
                        def parts = line.split(';')*.trim()
                        def slope = parts[0]
                        parts.drop(1).eachWithIndex { fmi_str, i ->
                            def fmi = fmi_str.toDouble()
                            if (fmi > best_fmi) {
                                best_fmi = fmi; best_slope = slope; best_intercept = intercepts[i]
                            }
                        }
                    }
                    [ best_slope, best_intercept ]
                }
            ch_fold_for_rffold = ch_fold_input
                .combine(ch_calibration)
                .map { meta, xmls, slope, intercept ->
                    [ meta + [ jackknife_slope: slope, jackknife_intercept: intercept ], xmls ]
                }
        } else {
            // Gate each fold group on all jackknife jobs for that cell_line completing.
            ch_fold_for_rffold = ch_fold_input
                .map { meta, xmls -> [ meta.id.toString(), meta, xmls ] }
                .join(
                    RNAFRAMEWORK_RFJACKKNIFE.out.csv
                        .map { meta, _csv -> [ meta.cell_line.toString(), 'done' ] }
                        .groupTuple()
                        .map { cell_line, _dones -> [ cell_line, 'done' ] }
                )
                .map { _id, meta, xmls, _done -> [ meta, xmls ] }
        }
    }

    RNAFRAMEWORK_RFFOLD (
        ch_fold_for_rffold
    )
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())

    //
    // Convert dot-plot outputs to .bp arc files for genome browser visualisation.
    // Genome-coordinate: remap via GTF.  Transcript-coordinate: use dot-plots directly.
    //
    def ch_dotplot_bp_input = RNAFRAMEWORK_RFFOLD.out.structures
        .combine(ch_reference_gtf_map)
        .flatMap { combined ->
            def meta      = combined[0]
            def fold_dir  = combined[1]
            def gtf_map   = combined[2]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def gtf_tuple = gtf_map[reference_key]
            if (!gtf_tuple) {
                log.warn("Skipping dotplot-to-bp conversion for '${reference_key}': no GTF available (expected for NCBI/viral references).")
                return []
            }
            return [ [ meta, fold_dir, gtf_tuple[1] ] ]
        }

    RNAFRAMEWORK_DOTPLOT2BP (
        ch_dotplot_bp_input,
        file("${projectDir}/bin/rnaframework_dotplot2bp.py", checkIfExists: true)
    )

    // Transcript-coordinate arc files — same dot-plots, skip GTF remapping.
    RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT (
        ch_dotplot_bp_input,
        file("${projectDir}/bin/rnaframework_dotplot2bp.py", checkIfExists: true)
    )

    MERGE_BP (
        RNAFRAMEWORK_DOTPLOT2BP.out.bp,
        file("${projectDir}/bin/merge_bp.py", checkIfExists: true)
    )

    MERGE_BP_TRANSCRIPT (
        RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT.out.bp,
        file("${projectDir}/bin/merge_bp.py", checkIfExists: true)
    )

    ch_versions = ch_versions.mix(RNAFRAMEWORK_DOTPLOT2BP.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BP.out.versions.first())

    emit:
    fold_input    = ch_fold_input                          // channel: [ val(meta), list(path(xml)) ] — grouped by cell_line
    structures    = RNAFRAMEWORK_RFFOLD.out.structures     // channel: [ val(meta), path(fold_dir) ]
    shannon_wig   = RNAFRAMEWORK_RFFOLD.out.shannon_wig    // channel: [ val(meta), path(wig) ]
    rffold_log    = RNAFRAMEWORK_RFFOLD.out.log            // channel: [ val(meta), path(log) ]
    jackknife_csv = ch_jackknife_csv                       // channel: [ val(meta), path(csv) ] — empty when jackknife not configured
    bp_dotplot    = RNAFRAMEWORK_DOTPLOT2BP.out.bp         // channel: [ val(meta), path(bp) ] — per-transcript, unmerged
    bp            = MERGE_BP.out.bp                        // channel: [ val(meta), path(bp) ] — merged per fold group
    bp_transcript = MERGE_BP_TRANSCRIPT.out.bp             // channel: [ val(meta), path(bp) ]
    versions      = ch_versions                            // channel: [ path(versions.yml) ]
}
