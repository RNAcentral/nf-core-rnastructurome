//
// FOLD_STRUCTURES — group rfnorm XML files by sample_group, optionally run rf-jackknife for
// slope/intercept calibration, run rf-fold to predict secondary structures, then convert
// dot-plot outputs to .bp arc files (genome- and transcript-coordinate) and merge them.
//

include { RNAFRAMEWORK_RFJACKKNIFE                                     } from '../../../modules/local/rnaframework/jackknife/main'
include { RNAFRAMEWORK_RFEVAL                                          } from '../../../modules/local/rnaframework/eval/main'
include { RNAFRAMEWORK_RFFOLD                                          } from '../../../modules/local/rnaframework/fold/main'
include { RNAFRAMEWORK_RFSTRUCTEXTRACT                                 } from '../../../modules/local/rnaframework/structextract/main'
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
    // Group rfnorm XMLs by sample_group so biological replicates are folded together.
    // The group key intentionally excludes replicate.
    //
    def ch_fold_input = ch_rfnorm_xml
        .map { meta, xml ->
            if (!meta.sample_group) {
                error("Missing sample_group for sample '${meta.id}'. rf-fold replicate grouping requires sample_group.")
            }
            def fold_group = meta.sample_group.toString()
            [ fold_group, [ meta, xml ] ]
        }
        .groupTuple()
        .map { fold_group, entries ->
            // Each entry is one rfnorm group (= one replicate) for this sample_group. rf-fold folds
            // replicates together via majority voting, taking one experiment DIRECTORY per replicate
            // (rf-fold exp1/ exp2/ ...), not a single merged dir. Sort by replicate for determinism,
            // then concatenate XMLs grouped by replicate so the module can rebuild per-replicate dirs
            // from the per-file input*/ staging using the aligned fold_replicate_sizes counts.
            def sortedEntries  = entries.sort { entry -> (entry[0].replicate ?: 'na').toString() }
            def base           = sortedEntries[0][0]
            def perRepXmls     = sortedEntries.collect { entry -> [ entry[1] ].flatten() }
            def orderedXmls    = perRepXmls.flatten()
            def replicateSizes = perRepXmls.collect { repXmls -> repXmls.size() }
            def replicates     = sortedEntries.collect { entry -> (entry[0].replicate ?: 'na').toString() }
            def sampleIds      = sortedEntries.collect { entry -> (entry[0].id ?: 'na').toString() }
            def foldMeta   = base + [
                id                    : fold_group,
                fold_group            : fold_group,
                fold_replicates       : replicates.join(','),
                fold_source_ids       : sampleIds.join(','),
                fold_experiment_count : sortedEntries.size(),   // replicate experiments; drives rf-fold -oc
                fold_replicate_sizes  : replicateSizes          // XMLs per replicate, aligned to orderedXmls
            ]
            [ foldMeta, orderedXmls ]
        }

    //
    // Optional rf-jackknife — quality assessment against a reference structure set.
    // When --jackknife_reference is provided, jackknife runs and rf-fold is gated on completion.
    // When --rfjackknife_pool_all is true, the optimal slope/intercept from the jackknife CSV
    // is parsed and injected directly into rf-fold meta, bypassing a separate calibration run.
    //
    def ch_fold_for_rffold = ch_fold_input
    def ch_jackknife_csv   = channel.empty()
    def ch_rfeval_csv      = channel.empty()

    if (pipeline_config.jackknife_reference) {
        def ch_jackknife_reference = channel.value(file(pipeline_config.jackknife_reference.toString(), checkIfExists: true))

        // Jackknife runs per rfnorm group (sample_group + replicate), not per fold group.
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
            // Gate each fold group on all jackknife jobs for that sample_group completing.
            ch_fold_for_rffold = ch_fold_input
                .map { meta, xmls -> [ meta.id.toString(), meta, xmls ] }
                .join(
                    RNAFRAMEWORK_RFJACKKNIFE.out.csv
                        .map { meta, _csv -> [ meta.sample_group.toString(), 'done' ] }
                        .groupTuple()
                        .map { sample_group, _dones -> [ sample_group, 'done' ] }
                )
                .map { _id, meta, xmls, _done -> [ meta, xmls ] }
        }
    }

    //
    // Optional rf-eval — evaluate agreement between reactivity data and a reference structure set.
    // Runs per rfnorm group (one XML set per sample_group+replicate) when --rfeval_reference is provided.
    //
    if (pipeline_config.rfeval_reference) {
        def ch_rfeval_reference = channel.value(file(pipeline_config.rfeval_reference.toString(), checkIfExists: true))

        RNAFRAMEWORK_RFEVAL (
            ch_rfnorm_xml,
            ch_rfeval_reference
        )
        ch_versions    = ch_versions.mix(RNAFRAMEWORK_RFEVAL.out.versions.first())
        ch_rfeval_csv  = RNAFRAMEWORK_RFEVAL.out.csv
    }

    def ch_fold_structures  = channel.empty()
    def ch_shannon_wig      = channel.empty()
    def ch_rffold_log       = channel.empty()
    def ch_bp_dotplot       = channel.empty()
    def ch_bp               = channel.empty()
    def ch_bp_transcript    = channel.empty()
    def ch_structextract    = channel.empty()

    if (!pipeline_config.stop_after_jackknife) {
        RNAFRAMEWORK_RFFOLD (
            ch_fold_for_rffold
        )
        ch_versions     = ch_versions.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())
        ch_fold_structures = RNAFRAMEWORK_RFFOLD.out.structures
        ch_shannon_wig     = RNAFRAMEWORK_RFFOLD.out.shannon_wig
        ch_rffold_log      = RNAFRAMEWORK_RFFOLD.out.log

        //
        // Optional rf-structextract — pull high-confidence, low-reactivity / low-Shannon structural
        // motifs out of the rf-fold output. Pairs each fold output directory (-ro) with the same
        // group's rf-norm XML reactivities (-xf). XML names repeat across replicates, so deduplicate
        // by filename (rf-structextract reads one reactivity profile per transcript).
        //
        if (pipeline_config.structextract) {
            def ch_structextract_input = RNAFRAMEWORK_RFFOLD.out.structures
                .map { meta, fold_dir -> [ meta.id.toString(), meta, fold_dir ] }
                .join(
                    ch_fold_for_rffold.map { meta, xmls -> [ meta.id.toString(), xmls ] }
                )
                .map { _id, meta, fold_dir, xmls ->
                    def seen    = [] as Set
                    def deduped = [xmls].flatten().findAll { x -> seen.add(x.name) }
                    [ meta, fold_dir, deduped ]
                }

            RNAFRAMEWORK_RFSTRUCTEXTRACT(ch_structextract_input)
            ch_versions      = ch_versions.mix(RNAFRAMEWORK_RFSTRUCTEXTRACT.out.versions.first())
            ch_structextract = RNAFRAMEWORK_RFSTRUCTEXTRACT.out.motifs
        }

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
        ch_versions = ch_versions.mix(RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT.out.versions.first())
        ch_versions = ch_versions.mix(MERGE_BP.out.versions.first())
        ch_versions = ch_versions.mix(MERGE_BP_TRANSCRIPT.out.versions.first())

        ch_bp_dotplot    = RNAFRAMEWORK_DOTPLOT2BP.out.bp
        ch_bp            = MERGE_BP.out.bp
        ch_bp_transcript = MERGE_BP_TRANSCRIPT.out.bp
    }

    emit:
    fold_input    = ch_fold_input       // channel: [ val(meta), list(path(xml)) ] — grouped by sample_group
    structures    = ch_fold_structures  // channel: [ val(meta), path(fold_dir) ]  — empty when stop_after_jackknife
    shannon_wig   = ch_shannon_wig      // channel: [ val(meta), path(wig) ]        — empty when stop_after_jackknife
    rffold_log    = ch_rffold_log       // channel: [ val(meta), path(log) ]        — empty when stop_after_jackknife
    jackknife_csv = ch_jackknife_csv    // channel: [ val(meta), path(csv) ]        — empty when --jackknife_reference not set
    rfeval_csv    = ch_rfeval_csv       // channel: [ val(meta), path(csv) ]        — empty when --rfeval_reference not set
    bp_dotplot    = ch_bp_dotplot       // channel: [ val(meta), path(bp) ]         — empty when stop_after_jackknife
    bp            = ch_bp              // channel: [ val(meta), path(bp) ]          — empty when stop_after_jackknife
    bp_transcript = ch_bp_transcript   // channel: [ val(meta), path(bp) ]          — empty when stop_after_jackknife
    structextract = ch_structextract   // channel: [ val(meta), path(dir) ]         — empty unless --structextract
    versions      = ch_versions                            // channel: [ path(versions.yml) ]
}
