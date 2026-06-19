/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { RNAFRAMEWORK_RFNORM      } from '../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFJACKKNIFE } from '../modules/local/rnaframework/jackknife/main'
include { RNAFRAMEWORK_RFFOLD      } from '../modules/local/rnaframework/fold/main'
include { RNAFRAMEWORK_DOTPLOT2BP                                  } from '../modules/local/dotplot2bp/main'
include { RNAFRAMEWORK_DOTPLOT2BP as RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT } from '../modules/local/dotplot2bp/main'
include { MERGE_BP                                                 } from '../modules/local/merge_bp/main'
include { MERGE_BP as MERGE_BP_TRANSCRIPT                          } from '../modules/local/merge_bp/main'
include { RNAFRAMEWORK_TORDAT    } from '../modules/local/tordat/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { FASTQ_QC_TRIM          } from '../subworkflows/local/fastq_qc_trim/main'
include { PREPARE_REFERENCES     } from '../subworkflows/local/prepare_references/main'
include { ALIGN_READS            } from '../subworkflows/local/align_reads/main'
include { QUANTIFY_REACTIVITY    } from '../subworkflows/local/quantify_reactivity/main'
include { VISUALISE_STRUCTURES   } from '../subworkflows/local/visualise/main'
include { BROWSER_TRACKS         } from '../subworkflows/local/browser_tracks/main'

// Pure helper functions (parsers, arg renderers, MultiQC table builders) — see rnastructurome_functions.nf
include {
    resolveReferenceKey
    resolveReferenceResolution
    uniqueReferenceResolution
    collectToMap
    buildStarAlignInputs
    parseFlagstatMappedReads
    parseRfcountCoveredTranscripts
    parseInferExperiment
    parseRfnormLog
    parseRffoldLog
    parseCutadaptCommandArg
    resolveRfNormNormMethod
    countProgressionMultiqc
    rfnormStatsMultiqc
    rffoldStatsMultiqc
    cutadaptAdaptersMultiqc
    filterSummaryParams
    addModuleOptionsSummary
} from './rnastructurome_functions.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASTRUCTUROME {

    take:
    ch_samplesheet         // channel: samplesheet read in from --input
    pipeline_config_input  // map: pipeline configuration captured at the entry workflow
    main:

    def pipeline_config = defaultPipelineConfig() + (pipeline_config_input ?: [:])
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // SUBWORKFLOW: FASTQ_QC_TRIM — cat → FastQC → UMI → cutadapt → FastQC
    //
    FASTQ_QC_TRIM (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM.out.multiqc_files)

    def ch_samplesheet_for_branching = FASTQ_QC_TRIM.out.reads_branched
    def ch_rtstop_trimmed_for_align = FASTQ_QC_TRIM.out.rtstop_trimmed
    def ch_map_trimmed_for_align    = FASTQ_QC_TRIM.out.map_trimmed

    //
    // SUBWORKFLOW: PREPARE_REFERENCES — resolve/fetch/sort/index reference artifacts
    //
    PREPARE_REFERENCES (
        ch_samplesheet_for_branching,
        pipeline_config
    )
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    def ch_reference_fasta_map          = PREPARE_REFERENCES.out.fasta_map
    def ch_reference_gtf_map            = PREPARE_REFERENCES.out.gtf_map
    def ch_all_reference_gtf            = PREPARE_REFERENCES.out.all_gtf
    def ch_reference_genome_fasta_keyed = PREPARE_REFERENCES.out.genome_fasta_keyed
    def ch_star_index                   = PREPARE_REFERENCES.out.star_index
    def ch_bowtie_index_map             = PREPARE_REFERENCES.out.bowtie_index_map
    def ch_bowtie2_index_map            = PREPARE_REFERENCES.out.bowtie2_index_map

    //
    // SUBWORKFLOW: ALIGN_READS — align → sort → dedup → stats → strandedness
    //
    ALIGN_READS (
        ch_rtstop_trimmed_for_align,
        ch_map_trimmed_for_align,
        ch_star_index,
        ch_all_reference_gtf,
        ch_bowtie_index_map,
        ch_bowtie2_index_map,
        ch_reference_fasta_map,
        pipeline_config
    )
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_READS.out.multiqc_files)

    def ch_markdup_bam_bai    = ALIGN_READS.out.markdup_bam_bai
    def ch_dedup_bam          = ALIGN_READS.out.dedup_bam
    def ch_strandedness_by_id = ALIGN_READS.out.strandedness_by_id

    //
    //
    // SUBWORKFLOW: QUANTIFY_REACTIVITY — rf-count(-genome) → rf-rctools → RC files
    //
    QUANTIFY_REACTIVITY (
        ch_markdup_bam_bai,
        ch_strandedness_by_id,
        ch_reference_genome_fasta_keyed,
        ch_reference_gtf_map,
        ch_reference_fasta_map,
        pipeline_config
    )
    ch_versions = ch_versions.mix(QUANTIFY_REACTIVITY.out.versions)

    def ch_rfcount_rc      = QUANTIFY_REACTIVITY.out.rc
    def ch_rfcount_rci     = QUANTIFY_REACTIVITY.out.rci
    def ch_rfcount_summary = QUANTIFY_REACTIVITY.out.summary
    def ch_rfcount_plots   = QUANTIFY_REACTIVITY.out.plots

    def ch_pre_dedup_mapped_reads = ALIGN_READS.out.flagstat_pre
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_post_dedup_mapped_reads = ALIGN_READS.out.flagstat_post
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_rfcount_covered_transcripts = ch_rfcount_summary
        .map { meta, summary_tsv -> [ meta.id.toString(), parseRfcountCoveredTranscripts(summary_tsv) ] }

    def ch_count_progression_mqc = ch_pre_dedup_mapped_reads
        .join(ch_post_dedup_mapped_reads)
        .join(ch_rfcount_covered_transcripts)
        .map { sample_id, mapped_pre, mapped_post, covered_transcripts ->
            def mappedPreLong = mapped_pre as long
            def mappedPostLong = mapped_post as long
            def removed = Math.max(mappedPreLong - mappedPostLong, 0L)
            def pctRemoved = mappedPreLong ? ((removed as double) / (mappedPreLong as double)) * 100.0d : 0.0d
            [
                sample_id,
                [
                    mapped_reads_pre_dedup     : mappedPreLong,
                    mapped_reads_post_dedup    : mappedPostLong,
                    pct_removed_by_dedup       : pctRemoved,
                    rfcount_covered_transcripts: covered_transcripts as long
                ]
            ]
        }
        .collect()
        .map { rows -> countProgressionMultiqc(rows) }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_count_progression_mqc.collectFile(
            name: 'count_progression_mqc.yaml',
            sort: true
        )
    )

    //
    // MODULE: rf-norm — normalise RC files to per-base reactivities (XML)
    //

    ch_rc_with_rci = ch_rfcount_rc
        .map { meta, rc -> [ meta.id.toString(), meta, rc ] }
        .join(
            ch_rfcount_rci
                .map { meta, rci -> [ meta.id.toString(), rci ] },
            remainder: true
        )
        .map { _sample_id, meta, rc, rci -> [ meta, rc, rci ?: [] ] }

    ch_rc_by_group = ch_rc_with_rci
        .map { meta, rc, rci ->
            if (!meta.cell_line || !meta.replicate) {
                error("Missing cell_line or replicate for sample '${meta.id}'. rf-norm requires both columns to pair samples safely.")
            }
            def group = "${meta.cell_line}_${meta.replicate}".toString()
            def condition = (meta.condition ?: 'treated').toLowerCase()
            [ group, condition, meta, rc, rci ]
        }

    ch_treated   = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'treated' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }
        .groupTuple()

    ch_untreated = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'untreated' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }

    ch_denatured = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'denatured' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }

    // Stage any available .rci sidecars alongside RC files so rf-norm can auto-discover them.
    ch_group_rci = ch_rc_by_group
        .filter  { _group, _condition, _meta, _rc, rci -> rci }
        .map     { group, _condition, _meta, _rc, rci -> [ group, rci ] }
        .groupTuple()

    // Enforce rf-norm pairing rules explicitly and select the representative meta per group.
    // Must be defined before the fallback channels so ch_treated_no_untreated can join against it.
    ch_group_meta = ch_rc_by_group
        .map     { group, condition, meta, _rc, _rci -> [ group, [ condition: condition, meta: meta ] ] }
        .groupTuple()
        .map     { group, entries ->
            def conditions = entries.collect { entry -> entry.condition }.toSet()
            if (!conditions.contains('treated') && conditions.contains('untreated')) {
                def offendingConditions = entries
                    .findAll { entry -> entry.condition == 'untreated' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }
                    .sort()
                    .join(', ')
                error("rf-norm requires a treated sample for cell_line/replicate group '${group}'. Invalid samples: ${offendingConditions}")
            }
            if (conditions.contains('denatured') && (!conditions.contains('treated') || !conditions.contains('untreated'))) {
                def offendingConditions = entries
                    .findAll { entry -> entry.condition == 'denatured' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }
                    .sort()
                    .join(', ')
                error("rf-norm requires treated and untreated samples for denatured controls in cell_line/replicate group '${group}'. Invalid samples: ${offendingConditions}")
            }
            [ group, entries[0].meta ]
        }

    // Fuzzy untreated pairing fallback via channel joins (avoids combine-with-empty-list issues).
    // When a treated group has no exact cell_line+replicate untreated match, the pipeline
    // falls back to an untreated sample that shares the same first cell_line base token
    // (e.g. MDA-MB-231_DMSO untreated covers MDA-MB-231_MTX treated) at the same replicate.
    //
    // Each untreated sample is emitted as [base_token, replicate, group, rc] for cross-matching.
    def ch_untreated_for_lookup = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'untreated' }
        .map     { group, _condition, meta, rc, _rci ->
            [ meta.cell_line.toString().tokenize('_')[0], meta.replicate.toString(), group, rc ]
        }

    // Treated groups with no direct untreated match — need a fallback.
    def ch_treated_no_untreated = ch_treated
        .join(ch_group_meta)
        .join(ch_untreated, remainder: true)
        .filter { _group, _treated_rcs, _base_meta, untreated_rc -> !untreated_rc }
        .map    { group, _treated_rcs, base_meta, _untreated_rc ->
            [ base_meta.cell_line.toString().tokenize('_')[0], base_meta.replicate.toString(), group ]
        }

    // Cross-product treated-without-untreated × available untreated, filter on base+rep match,
    // then group by treated group to validate uniqueness before selecting the fallback.
    // When ch_untreated_for_lookup is empty (all-treated datasets), this channel is empty too
    // and the remainder:true join below correctly produces null for untreated.
    def ch_fallback_untreated = ch_treated_no_untreated
        .combine(ch_untreated_for_lookup)
        .filter { treated_base, treated_rep, _group, unt_base, unt_rep, _unt_group, _unt_rc ->
            treated_base == unt_base && treated_rep == unt_rep
        }
        .map { _treated_base, _treated_rep, group, _unt_base, _unt_rep, unt_group, unt_rc ->
            [ group, [ untreated_group: unt_group, rc: unt_rc ] ]
        }
        .groupTuple()
        .map { group, candidates ->
            if (candidates.size() > 1) {
                def candidateGroups = candidates.collect { c -> c.untreated_group }.sort().join(', ')
                error("Ambiguous untreated fallback for '${group}': multiple untreated groups share the same cell_line base and replicate: ${candidateGroups}.")
            }
            log.warn "No exact untreated match for '${group}' — falling back to '${candidates[0].untreated_group}' (shared cell_line base at same replicate)."
            [ group, candidates[0].rc ]
        }

    // Resolved untreated: direct exact match OR fuzzy fallback.
    def ch_resolved_untreated = ch_untreated.mix(ch_fallback_untreated)

    ch_norm_input = ch_treated
        .join(ch_group_meta)
        .join(ch_resolved_untreated, remainder: true)
        .join(ch_denatured, remainder: true)
        .join(ch_group_rci, remainder: true)
        .map { group, treated_rcs, base_meta, untreated_rc, denatured_rc, rci_files ->
            def hasUntreated = untreated_rc ? true : false
            def hasDenatured = denatured_rc ? true : false
            def principle    = (base_meta.principle ?: '').toLowerCase()
            def scoringMethod = principle == 'map'
                ? (hasUntreated ? 3 : 4)
                : (hasUntreated ? 1 : 2)
            def normMethod = resolveRfNormNormMethod(pipeline_config, scoringMethod)
            def gmeta = base_meta + [
                id                    : group,
                rfnorm_has_untreated  : hasUntreated,
                rfnorm_has_denatured  : hasDenatured,
                rfnorm_scoring_method : scoringMethod,
                rfnorm_norm_method    : normMethod
            ]
            [ gmeta, treated_rcs, untreated_rc ?: [], denatured_rc ?: [], rci_files ?: [] ]
        }

    RNAFRAMEWORK_RFNORM (
        ch_norm_input
    )

    //
    // MODULE: rf-fold — predict RNA secondary structures from reactivity XML
    //
    // Fold across all available replicate XMLs per experimental group.
    // Group key intentionally excludes replicate so biological replicates can
    // be folded together when present.
    ch_fold_input = RNAFRAMEWORK_RFNORM.out.xml
        .map { meta, xml ->
            if (!meta.cell_line) {
                error("Missing cell_line for sample '${meta.id}'. rf-fold replicate grouping requires cell_line.")
            }
            def fold_group = meta.cell_line.toString()
            [ fold_group, [ meta, xml ] ]
        }
        .groupTuple()
        .map { fold_group, entries ->
            def metas = entries.collect { entry -> entry[0] }
            def xmls = entries.collect { entry -> entry[1] }.flatten()
            def base = metas[0]
            def replicates = metas.collect { meta -> (meta.replicate ?: 'na').toString() }.unique().sort()
            def sampleIds = metas.collect { meta -> (meta.id ?: 'na').toString() }.unique().sort()
            def foldMeta = base + [
                id                 : fold_group,
                fold_group         : fold_group,
                fold_replicates    : replicates.join(','),
                fold_source_ids    : sampleIds.join(','),
                fold_xml_count     : xmls.size()
            ]
            [ foldMeta, xmls ]
        }

    //
    // MODULE: rf-jackknife — optional normalisation quality assessment against a reference structure set.
    // When --jackknife_reference is provided, rf-jackknife runs and rf-fold is gated on its completion.
    // When --rfjackknife_pool_all is true, the optimal slope/intercept from the jackknife CSV is parsed
    // and injected directly into rf-fold (via meta), bypassing the need for a separate calibration run.
    // When omitted, rf-fold runs directly using --rffold_slope/--rffold_intercept if set.
    //
    def ch_fold_for_rffold = ch_fold_input

    if (pipeline_config.jackknife_reference) {
        def ch_jackknife_reference = channel.value(file(pipeline_config.jackknife_reference.toString(), checkIfExists: true))

        // Jackknife runs per rfnorm group (cell_line + replicate), not per fold group.
        // Fold groups flatten XMLs from multiple replicates which produce identically-named
        // files (e.g. 16S_rRNA.xml); using input*/* staging gives each file its own
        // numbered directory so rf-jackknife receives them as separate experiment dirs.
        def ch_jackknife_input
        if (pipeline_config.rfjackknife_pool_all as Boolean) {
            ch_jackknife_input = RNAFRAMEWORK_RFNORM.out.xml
                .flatMap { _meta, xmls -> [xmls].flatten() }
                .collect()
                .map { all_xmls -> [ [ id: 'all_groups', fold_group: 'all_groups' ], all_xmls ] }
        } else {
            ch_jackknife_input = RNAFRAMEWORK_RFNORM.out.xml
        }

        RNAFRAMEWORK_RFJACKKNIFE (
            ch_jackknife_input,
            ch_jackknife_reference
        )
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFJACKKNIFE.out.versions.first())

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

    def ch_dotplot_bp_input = RNAFRAMEWORK_RFFOLD.out.structures
        .combine(ch_reference_gtf_map)
        .flatMap { combined ->
            def meta = combined[0]
            def fold_dir = combined[1]
            def gtf_map = combined[2]
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

    //
    // MODULE: merge_bp — merge per-transcript .bp files into a single file per fold group for genome browser visualisation
    //
    MERGE_BP (
        RNAFRAMEWORK_DOTPLOT2BP.out.bp,
        file("${projectDir}/bin/merge_bp.py", checkIfExists: true)
    )

    MERGE_BP_TRANSCRIPT (
        RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT.out.bp,
        file("${projectDir}/bin/merge_bp.py", checkIfExists: true)
    )

    //
    // SUBWORKFLOW: VISUALISE_STRUCTURES — R2DT / ViennaRNA 2D structure diagrams
    //
    VISUALISE_STRUCTURES (
        RNAFRAMEWORK_RFNORM.out.xml,
        RNAFRAMEWORK_RFFOLD.out.structures,
        ch_fold_input,
        ch_reference_fasta_map,
        pipeline_config
    )
    ch_versions = ch_versions.mix(VISUALISE_STRUCTURES.out.versions)

    //
    // SUBWORKFLOW: BROWSER_TRACKS — rf-wiggle → BigWig genome/transcript tracks (+ Shannon)
    //
    BROWSER_TRACKS (
        RNAFRAMEWORK_RFNORM.out.xml,
        RNAFRAMEWORK_RFFOLD.out.shannon_wig,
        ch_reference_gtf_map,
        pipeline_config
    )
    ch_versions = ch_versions.mix(BROWSER_TRACKS.out.versions)

    //
    // MODULE: tordat — compile rf-norm XML + rf-fold .db structures into RDAT format
    //
    // Build per-reference source file name channels so the RDAT COMMENT records
    // the exact Ensembl filenames (e.g. Homo_sapiens.GRCh38.114.cdna.all.fa.gz)
    // or NCBI accessions (e.g. EU081230.1) rather than the generic pipeline names.
    def ch_ref_fasta_names
    def ch_ref_gtf_names
    if (pipeline_config.fasta) {
        def local_fasta_name = file(pipeline_config.fasta.toString()).name
        def local_gtf_name   = pipeline_config.gtf ? file(pipeline_config.gtf.toString()).name : ''
        ch_ref_fasta_names = PREPARE_REFERENCES.out.local_fasta_sorted
            .map { meta, _fasta -> [ meta.id.toString(), local_fasta_name ] }
        ch_ref_gtf_names = PREPARE_REFERENCES.out.local_fasta_sorted
            .map { meta, _fasta -> [ meta.id.toString(), local_gtf_name ] }
    } else {
        ch_ref_fasta_names = PREPARE_REFERENCES.out.ensembl_fasta_source_url
            .map { meta, urls_file ->
                def names = urls_file.readLines().findAll { line -> line.trim() }
                    .collect { line -> line.tokenize('/').last() }.join(' + ')
                [ meta.id.toString(), names ]
            }
            .mix(PREPARE_REFERENCES.out.ncbi_source_accessions
                .map { meta, acc_file ->
                    def accs = acc_file.readLines()
                        .findAll { line -> line.trim() && !line.startsWith('stub:') }
                        .collect { line -> line.tokenize('/').last() }.join(', ')
                    [ meta.id.toString(), accs ?: meta.id.toString() ]
                })
        ch_ref_gtf_names = PREPARE_REFERENCES.out.ensembl_gtf_source_urls
            .map { meta, urls_file ->
                def name = urls_file.readLines().find { line -> line.trim() }?.tokenize('/')?.last() ?: ''
                [ meta.id.toString(), name ]
            }
            .mix(PREPARE_REFERENCES.out.gtf_local
                .map { meta, gtf_file -> [ meta.id.toString(), gtf_file.name ] })
            .mix(PREPARE_REFERENCES.out.ncbi_source_accessions
                .map { meta, _acc -> [ meta.id.toString(), '' ] })
    }

    def ch_rdat_input = ch_fold_input
        .map { meta, xml -> [ meta.id.toString(), meta, xml ] }
        .combine(
            RNAFRAMEWORK_RFFOLD.out.structures
                .map { meta, fold_dir -> [ meta.id.toString(), fold_dir ] },
            by: 0
        )
        .map { _key, fold_meta, xml, fold_dir ->
            def reference_key = resolveReferenceKey(fold_meta, pipeline_config.organism)
            [ reference_key.toString(), fold_meta, xml, fold_dir ]
        }
        .combine(ch_ref_fasta_names, by: 0)
        .combine(ch_ref_gtf_names, by: 0)
        .map { _ref_key, fold_meta, xml, fold_dir, fasta_name, gtf_name ->
            [ fold_meta + [ fasta_name: fasta_name, gtf_name: gtf_name ], xml, fold_dir ]
        }

    RNAFRAMEWORK_TORDAT (
        ch_rdat_input,
        file("${projectDir}/bin/rnaframework_to_rdat.py", checkIfExists: true)
    )

    // Add RNAframework outputs to MultiQC input collection.
    ch_multiqc_files = ch_multiqc_files.mix(ch_rfcount_rc.collect { rc_file -> rc_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_rfcount_plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.xml.collect { xml_file -> xml_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFFOLD.out.structures.collect { fold_dir -> fold_dir[1] })

    // RF-norm summary table: one row per normalisation group (cell_line + replicate).
    // Join RF-count covered transcript counts per group (max across treated replicates).
    def ch_rfcount_covered_by_group = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'treated' }
        .map     { group, _condition, meta, _rc, _rci -> [ meta.id.toString(), group ] }
        .join(ch_rfcount_covered_transcripts)
        .map     { _sample_id, group, covered -> [ group, covered as long ] }
        .groupTuple()
        .map     { group, covered_list -> [ group, covered_list.max() ] }

    def ch_rfnorm_stats_mqc = RNAFRAMEWORK_RFNORM.out.log
        .map { meta, log -> [ meta.id.toString(), parseRfnormLog(log) ] }
        .join(ch_rfcount_covered_by_group, remainder: true)
        .map { group_id, rfnorm_stats, rfcount_covered ->
            [ group_id, [ rfcount_covered: (rfcount_covered ?: 0L) as long, covered: rfnorm_stats?.covered ?: 0L ] ]
        }
        .collect()
        .map { rows -> rfnormStatsMultiqc(rows) }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_rfnorm_stats_mqc.collectFile(name: 'rfnorm_stats_mqc.yaml', sort: true)
    )

    // RF-fold summary table: one row per fold group (cell_line; may span replicates).
    def ch_rffold_stats_mqc = RNAFRAMEWORK_RFFOLD.out.log
        .map { meta, log -> [ meta.id.toString(), parseRffoldLog(log) ] }
        .collect()
        .map { rows -> rffoldStatsMultiqc(rows) }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_rffold_stats_mqc.collectFile(name: 'rffold_stats_mqc.yaml', sort: true)
    )

    // RNAframework module versions collected via ch_versions below

    //
    // MODULE: multiqc — aggregate pipeline quality control reports
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = pipeline_config.multiqc_config ?
        channel.fromPath(pipeline_config.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = pipeline_config.multiqc_logo ?
        channel.fromPath(pipeline_config.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    summary_params      = filterSummaryParams(summary_params)
    summary_params      = addModuleOptionsSummary(summary_params, pipeline_config)
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = pipeline_config.multiqc_methods_description ?
        file(pipeline_config.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description, pipeline_config))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect().map { f -> [f] }
            .combine(ch_multiqc_config.mix(ch_multiqc_custom_config).collect().map { c -> [c] })
            .combine(ch_multiqc_logo.collect().ifEmpty([]).map { l -> [l] })
            .map { files, config, logo -> [ [ id: 'multiqc' ], files, config, logo, [], [] ] }
    )

    //
    // Collect software versions from all modules.
    // Modules using `path "versions.yml"` (old pattern) must be mixed in explicitly.
    // Modules using `topic: versions` (new pattern) are collected automatically by
    // channel.topic("versions") below and do NOT need to be listed here.
    //   topic-pattern: CUTADAPT_*, BOWTIE2_*, UMITOOLS_*, all SAMTOOLS_* modules
    //   file-pattern:  FASTQC, BOWTIE_BUILD/ALIGN, local modules
    //
    // FASTQC, SAMtools, cutadapt, bowtie2, umitools use topic: versions → captured by channel.topic("versions") below
    // Old-style modules (emit: versions) must be mixed in explicitly
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_DOTPLOT2BP.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_DOTPLOT2BP_TRANSCRIPT.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BP.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BP_TRANSCRIPT.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_TORDAT.out.versions.first())

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${pipeline_config.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnastructurome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    multiqc_report   = MULTIQC.out.report.toList()             // channel: /path/to/multiqc_report.html
    mapped_bam       = ch_dedup_bam                            // channel: [ val(meta), path(bam) ]
    normalized_xml   = RNAFRAMEWORK_RFNORM.out.xml             // channel: [ val(meta), path(xml) ]
    jackknife_csv    = pipeline_config.jackknife_reference ? RNAFRAMEWORK_RFJACKKNIFE.out.csv : channel.empty()  // channel: [ val(meta), path(csv) ]
    fold_structures  = RNAFRAMEWORK_RFFOLD.out.structures      // channel: [ val(meta), path(dir) ]
    fold_bp          = RNAFRAMEWORK_DOTPLOT2BP.out.bp          // channel: [ val(meta), path(bp) ]
    merged_bp        = MERGE_BP.out.bp                        // channel: [ val(meta), path(*_merged.bp) ]
    versions         = ch_versions                             // channel: [ path(versions.yml) ]

}


def defaultPipelineConfig() {
    [
        organism                          : null,
        fasta                             : null,
        genome_fasta                      : null,
        transcriptome_fasta               : null,
        gtf                               : null,
        transcriptome                     : false,
        genomes                           : null,
        ensembl_species_map               : [
            'human': 'homo_sapiens',
            'mouse': 'mus_musculus',
            'yeast': 'saccharomyces_cerevisiae',
            'rat'  : 'rattus_norvegicus'
        ],
        ensembl_release                   : 'current',
        ensembl_base_url                  : 'https://ftp.ensembl.org/pub',
        ncbi_accessions_map               : [:],
        multiqc_config                    : null,
        multiqc_logo                      : null,
        multiqc_methods_description       : null,
        outdir                            : null,
        input                             : null,
        umi_pattern                       : null,
        bowtie_manual_only                : false,
        bowtie_mapping_params             : null,
        bowtie_k                          : null,
        bowtie_all                        : false,
        bowtie_norc                       : false,
        bowtie_trim5                      : 0,
        bowtie_trim3                      : 0,
        bowtie_seedlen                    : null,
        bowtie_n                          : 2,
        bowtie_v                          : null,
        bowtie_max                        : 1,
        bowtie_chunkmbs                   : 128,
        bowtie2_N                         : 0,
        bowtie2_D                         : 15,
        bowtie2_R                         : 2,
        bowtie2_mp                        : '6,2',
        bowtie2_dpad                      : 15,
        bowtie2_rdg                       : '5,3',
        bowtie2_rfg                       : '5,3',
        bowtie2_softclip                  : false,
        bowtie2_ma                        : 2,
        bowtie2_dovetail                  : false,
        jackknife_reference               : null,
        rfjackknife_pool_all              : true,
        rfnorm_reactive_bases             : null,
        rfnorm_remap_reactivities         : false,
        rfnorm_norm_window                : null,
        rfnorm_window_offset              : null,
        rfnorm_dynamic_window             : null,
        rfnorm_norm_independent           : false,
        rfnorm_norm_factor                : null,
        rfnorm_norm_method                : null,
        rfnorm_raw                        : false,
        rfnorm_pseudocount                : null,
        rfnorm_max_score                  : null,
        rfnorm_ignore_lower_than_untreated: false,
        rfnorm_max_untreated_mut          : null,
        rfnorm_max_mutation_rate          : null,
        rfnorm_mean_coverage              : 0,
        rfnorm_median_coverage            : 0,
        rfnorm_nan                        : 1000,
        rnaframework_r_path               : '/usr/bin/R'
    ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
