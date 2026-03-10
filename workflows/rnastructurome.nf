/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_PRE   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST  } from '../modules/nf-core/fastqc/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { CUTADAPT as CUTADAPT_RTSTOP } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_MAP    } from '../modules/nf-core/cutadapt/main'
include { UMITOOLS_EXTRACT       } from '../modules/nf-core/umitools/extract/main'
include { UMITOOLS_DEDUP         } from '../modules/nf-core/umitools/dedup/main'
include { BOWTIE_BUILD          } from '../modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN          } from '../modules/nf-core/bowtie/align/main'
include { BOWTIE2_BUILD         } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN         } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT         } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SORT    } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL   } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MARKDUP      } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_STATS        } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT     } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_PRE } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS     } from '../modules/nf-core/samtools/idxstats/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_FAIDX        } from '../modules/nf-core/samtools/faidx/main'
include { RNAFRAMEWORK_RFCOUNT  } from '../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFNORM   } from '../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFFOLD   } from '../modules/local/rnaframework/fold/main'
include { ENSEMBL_TRANSCRIPTOME } from '../modules/local/ensembl/transcriptome/main'
include { FASTA_SORT            } from '../modules/local/fasta/sort/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnastructurome_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASTRUCTUROME {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    def ch_samplesheet_checked = ch_samplesheet.map { meta, reads ->
        def principle = (meta.principle ?: '').toLowerCase()
        if (!(principle in ['rt-stop', 'map'])) {
            error("Unsupported principle '${meta.principle}' for sample '${meta.id}'. Expected one of: RT-stop, MaP.")
        }
        [meta, reads]
    }

    //
    // MODULE: cat/fastq — merge resequenced FASTQ files per sample before QC
    //
    CAT_FASTQ (
        ch_samplesheet_checked
    )

    def ch_pretrim_fastqc_input = CAT_FASTQ.out.reads
    def ch_samplesheet_for_branching = CAT_FASTQ.out.reads

    //
    // MODULE: fastqc (pre-trim) — quality control on raw reads
    //
    FASTQC_PRE (
        ch_pretrim_fastqc_input
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    // Branch by probing principle so RT-stop and MaP can use different default cutadapt args
    def principle_branches = ch_samplesheet_for_branching.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    def ch_reads_for_umi = principle_branches.rtstop.mix(principle_branches.map)
    def ch_reads_without_umi = ch_reads_for_umi.filter { meta, _reads ->
        !((meta.umi_pattern ?: '').toString().trim())
    }

    //
    // MODULE: umi_tools extract — extract UMIs from reads
    //
    UMITOOLS_EXTRACT (
        ch_reads_for_umi
    )

    def ch_reads_after_umi = ch_reads_without_umi.mix(UMITOOLS_EXTRACT.out.reads)
    def reads_for_cutadapt = ch_reads_after_umi.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }
    def ch_rtstop_reads_for_cutadapt = reads_for_cutadapt.rtstop
    def ch_map_reads_for_cutadapt = reads_for_cutadapt.map
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.log.collect { umi_log -> umi_log[1] })

    //
    // MODULE: cutadapt (RT-stop) — trim RT-stop reads
    //
    CUTADAPT_RTSTOP (
        ch_rtstop_reads_for_cutadapt
    )

    //
    // MODULE: cutadapt (MaP) — trim MaP reads
    //
    CUTADAPT_MAP (
        ch_map_reads_for_cutadapt
    )

    ch_trimmed_reads = CUTADAPT_RTSTOP.out.reads.mix(CUTADAPT_MAP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_RTSTOP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_MAP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    def ch_cutadapt_adapter_mqc = CUTADAPT_RTSTOP.out.log
        .map { meta, cutadapt_log ->
            [
                meta.id.toString(),
                [
                    cutadapt_mode: 'RT-stop',
                    adapter_5p  : parseCutadaptCommandArg(cutadapt_log, '-g'),
                    adapter_3p  : parseCutadaptCommandArg(cutadapt_log, '-a')
                ]
            ]
        }
        .mix(
            CUTADAPT_MAP.out.log.map { meta, cutadapt_log ->
                [
                    meta.id.toString(),
                    [
                        cutadapt_mode: 'MaP',
                        adapter_5p  : parseCutadaptCommandArg(cutadapt_log, '-g'),
                        adapter_3p  : parseCutadaptCommandArg(cutadapt_log, '-a')
                    ]
                ]
            }
        )
        .collect()
        .map { rows -> cutadaptAdaptersMultiqc(rows) }
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_cutadapt_adapter_mqc.collectFile(
            name: 'cutadapt_adapters_mqc.yaml',
            sort: true
        )
    )

    //
    // MODULE: fastqc (post-trim) — quality control on trimmed reads
    //
    FASTQC_POST (
        ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    ch_reference_requests = ch_samplesheet_for_branching
        .map { meta, _reads ->
            def reference_key = resolveReferenceKey(meta, params.organism)
            if (params.fasta) {
                return [ reference_key, "path::${params.fasta.toString()}" ]
            }

            def genome_entry = params.genomes?.containsKey(reference_key) ? params.genomes[reference_key] : null
            def transcript_fasta = genome_entry?.transcript_fasta ?: genome_entry?.transcriptome ?: genome_entry?.cdna
            if (transcript_fasta) {
                return [ reference_key, "path::${transcript_fasta.toString()}" ]
            }

            def ensembl_species = genome_entry?.ensembl_species ?: params.ensembl_species_map?.get(reference_key) ?: (reference_key ==~ /[a-z]+_[a-z0-9_]+/ ? reference_key : null)
            if (!ensembl_species) {
                error("No transcript FASTA resolved for reference '${reference_key}'. Provide --fasta with a transcript FASTA, set params.genomes['${reference_key}'].transcript_fasta (or transcriptome/cdna), or set params.genomes['${reference_key}'].ensembl_species.")
            }
            [ reference_key, "ensembl::${ensembl_species.toLowerCase()}" ]
        }
        .groupTuple()
        .map { reference_key, resolutions ->
            def unique_resolutions = resolutions.unique()
            if (unique_resolutions.size() != 1) {
                error("Multiple transcript reference resolutions were detected for reference '${reference_key}': ${unique_resolutions.join(', ')}")
            }
            [ reference_key, unique_resolutions[0] ]
        }

    ch_reference_local = ch_reference_requests
        .filter { _reference_key, resolution -> resolution.startsWith('path::') }
        .map { reference_key, resolution ->
            def fasta_path = resolution - 'path::'
            [ [ id: reference_key, organism: reference_key ], file(fasta_path, checkIfExists: true) ]
        }

    ch_reference_ensembl_input = ch_reference_requests
        .filter { _reference_key, resolution -> resolution.startsWith('ensembl::') }
        .map { reference_key, resolution ->
            def ensembl_species = resolution - 'ensembl::'
            [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species ], ensembl_species ]
        }

    ENSEMBL_TRANSCRIPTOME (
        ch_reference_ensembl_input
    )
    ch_versions = ch_versions.mix(ENSEMBL_TRANSCRIPTOME.out.versions)

    ch_all_reference_fasta = ch_reference_local.mix(ENSEMBL_TRANSCRIPTOME.out.fasta)

    FASTA_SORT (
        ch_all_reference_fasta
    )
    ch_versions = ch_versions.mix(FASTA_SORT.out.versions)

    ch_reference_fasta_keyed = FASTA_SORT.out.fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }

    ch_rtstop_reference_fasta = principle_branches.rtstop
        .map { meta, _reads -> [ resolveReferenceKey(meta, params.organism), true ] }
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, _flag, ref_tuple -> ref_tuple }

    ch_map_reference_fasta = principle_branches.map
        .map { meta, _reads -> [ resolveReferenceKey(meta, params.organism), true ] }
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, _flag, ref_tuple -> ref_tuple }

    //
    // MODULE: bowtie-build — build Bowtie v1 indices for RT-stop alignment
    //
    BOWTIE_BUILD (
        ch_rtstop_reference_fasta
    )

    //
    // MODULE: bowtie2-build — build Bowtie2 indices for MaP alignment
    //
    BOWTIE2_BUILD (
        ch_map_reference_fasta
    )

    ch_bowtie_index_keyed = BOWTIE_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
    ch_bowtie2_index_keyed = BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }

    ch_rtstop_align_inputs = CUTADAPT_RTSTOP.out.reads
        .map { meta, reads ->
            [ resolveReferenceKey(meta, params.organism), [meta, reads] ]
        }
        .join(ch_bowtie_index_keyed)
        .map { _reference_key, reads_tuple, index_tuple -> [ reads_tuple, index_tuple ] }

    ch_map_align_inputs = CUTADAPT_MAP.out.reads
        .map { meta, reads ->
            [ resolveReferenceKey(meta, params.organism), [meta, reads] ]
        }
        .join(ch_bowtie2_index_keyed)
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, reads_tuple, index_tuple, fasta_tuple -> [ reads_tuple, index_tuple, fasta_tuple ] }

    //
    // MODULE: bowtie align — align RT-stop reads with Bowtie v1
    //
    BOWTIE_ALIGN (
        ch_rtstop_align_inputs.map { reads_tuple, _index_tuple -> reads_tuple },
        ch_rtstop_align_inputs.map { _reads_tuple, index_tuple -> index_tuple },
        false
    )

    //
    // MODULE: bowtie2 align — align MaP reads with Bowtie2
    //
    BOWTIE2_ALIGN (
        ch_map_align_inputs.map { reads_tuple, _index_tuple, _fasta_tuple -> reads_tuple },
        ch_map_align_inputs.map { _reads_tuple, index_tuple, _fasta_tuple -> index_tuple },
        ch_map_align_inputs.map { _reads_tuple, _index_tuple, fasta_tuple -> fasta_tuple },
        false,
        false
    )

    ch_mapped_bam = BOWTIE_ALIGN.out.bam.mix(BOWTIE2_ALIGN.out.bam)
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect { bowtie_log -> bowtie_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { bowtie2_log -> bowtie2_log[1] })

    ch_sorted_inputs = ch_mapped_bam
        .map { meta, bam ->
            [ resolveReferenceKey(meta, params.organism), [meta, bam] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, bam_tuple, fasta_tuple ->
            [ bam_tuple, fasta_tuple ]
        }

    //
    // MODULE: samtools sort — coordinate-sort mapped BAMs
    //
    SAMTOOLS_SORT (
        ch_sorted_inputs.map { bam_tuple, _fasta_tuple -> bam_tuple },
        ch_sorted_inputs.map { _bam_tuple, fasta_tuple -> fasta_tuple },
        false
    )

    //
    // MODULE: samtools index (sorted) — index sorted BAMs
    //
    SAMTOOLS_INDEX_SORT (
        SAMTOOLS_SORT.out.bam
    )

    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_SORT.out.bai.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
        .map { _sample_id, bam_tuple, bai_tuple ->
            [ bam_tuple[0], bam_tuple[1], bai_tuple[1] ]
        }

    def dedup_branches = ch_sorted_bam_bai.branch { meta, _bam, _bai ->
        umi:     (meta.umi_pattern ?: '').toString().trim()
        non_umi: !((meta.umi_pattern ?: '').toString().trim())
    }

    //
    // MODULE: samtools flagstat — collect pre-dedup flag statistics
    //
    SAMTOOLS_FLAGSTAT_PRE (
        ch_sorted_bam_bai
    )

    //
    // MODULE: umi_tools dedup — deduplicate UMI-tagged BAMs
    //
    UMITOOLS_DEDUP (
        dedup_branches.umi,
        false
    )

    //
    // MODULE: samtools faidx — index reference FASTA for markdup
    //
    SAMTOOLS_FAIDX (
        FASTA_SORT.out.fasta.map { meta, fasta -> [meta, fasta, []] },
        false
    )

    ch_reference_fasta_fai_keyed = SAMTOOLS_FAIDX.out.fai
        .map { meta, fai -> [ meta.id.toString(), [meta, fai] ] }
        .join(ch_reference_fasta_keyed)
        .map { reference_key, fai_tuple, fasta_tuple ->
            [ reference_key, [fasta_tuple[0], fasta_tuple[1], fai_tuple[1]] ]
        }

    ch_markdup_inputs = dedup_branches.non_umi
        .map { meta, bam, _bai ->
            [ resolveReferenceKey(meta, params.organism), [meta, bam] ]
        }
        .join(ch_reference_fasta_fai_keyed)
        .map { _reference_key, bam_tuple, fasta_fai_tuple ->
            [ bam_tuple, fasta_fai_tuple ]
        }

    //
    // MODULE: samtools markdup — deduplicate non-UMI BAMs
    //
    SAMTOOLS_MARKDUP (
        ch_markdup_inputs.map { bam_tuple, _fasta_fai_tuple -> bam_tuple },
        ch_markdup_inputs.map { _bam_tuple, fasta_fai_tuple -> fasta_fai_tuple }
    )

    ch_dedup_bam = UMITOOLS_DEDUP.out.bam.mix(SAMTOOLS_MARKDUP.out.bam)

    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_DEDUP.out.log.collect { dedup_log -> dedup_log[1] })

    //
    // MODULE: samtools index (final) — index deduplicated BAMs
    //
    SAMTOOLS_INDEX_FINAL (
        ch_dedup_bam
    )

    ch_markdup_bam_bai = ch_dedup_bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_FINAL.out.bai.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
        .map { _sample_id, bam_tuple, bai_tuple ->
            [ bam_tuple[0], bam_tuple[1], bai_tuple[1] ]
        }

    ch_stats_inputs = ch_markdup_bam_bai
        .map { meta, bam, bai ->
            [ resolveReferenceKey(meta, params.organism), [meta, bam, bai] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, bam_bai_tuple, fasta_tuple ->
            [ bam_bai_tuple, fasta_tuple ]
        }

    //
    // MODULE: samtools stats — collect alignment statistics
    //
    SAMTOOLS_STATS (
        ch_stats_inputs.map { bam_bai_tuple, _fasta_tuple -> bam_bai_tuple },
        ch_stats_inputs.map { _bam_bai_tuple, fasta_tuple -> fasta_tuple }
    )

    //
    // MODULE: samtools flagstat — collect flag statistics
    //
    SAMTOOLS_FLAGSTAT (
        ch_markdup_bam_bai
    )

    //
    // MODULE: samtools idxstats — collect per-reference mapping statistics
    //
    SAMTOOLS_IDXSTATS (
        ch_markdup_bam_bai
    )

    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { stats_file -> stats_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect { flagstat_file -> flagstat_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect { idxstats_file -> idxstats_file[1] })

    //
    // Collect software versions for MultiQC.
    //
    ch_versions_for_multiqc_files = channel.empty()
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(FASTQC_PRE.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(FASTQC_POST.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(BOWTIE_BUILD.out.versions)
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(BOWTIE_ALIGN.out.versions)

    ch_versions_for_multiqc_tuples = channel.empty()
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(CAT_FASTQ.out.versions_cat)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(UMITOOLS_EXTRACT.out.versions_umitools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(UMITOOLS_DEDUP.out.versions_umitools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(CUTADAPT_RTSTOP.out.versions_cutadapt)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(CUTADAPT_MAP.out.versions_cutadapt)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_BUILD.out.versions_bowtie2)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_bowtie2)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_pigz)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_SORT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_INDEX_SORT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_MARKDUP.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_INDEX_FINAL.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_STATS.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_FLAGSTAT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_IDXSTATS.out.versions_samtools)

    //
    // MODULE: rf-count — per-base RT-stop or mutation counts from deduplicated BAM
    //
    ch_rfcount_with_fasta = ch_markdup_bam_bai
        .map { meta, bam, bai ->
            [ resolveReferenceKey(meta, params.organism), [meta, bam, bai] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _reference_key, bam_bai_tuple, fasta_tuple ->
            [ bam_bai_tuple, fasta_tuple ]
        }

    RNAFRAMEWORK_RFCOUNT (
        ch_rfcount_with_fasta.map { bam_bai_tuple, _fasta_tuple -> bam_bai_tuple },
        ch_rfcount_with_fasta.map { _bam_bai_tuple, fasta_tuple -> fasta_tuple }
    )

    def ch_pre_dedup_mapped_reads = SAMTOOLS_FLAGSTAT_PRE.out.flagstat
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_post_dedup_mapped_reads = SAMTOOLS_FLAGSTAT.out.flagstat
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_rfcount_covered_transcripts = RNAFRAMEWORK_RFCOUNT.out.summary
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
                    reads_removed_by_dedup     : removed,
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

    ch_rc_with_rci = RNAFRAMEWORK_RFCOUNT.out.rc
        .map { meta, rc -> [ meta.id.toString(), meta, rc ] }
        .join(
            RNAFRAMEWORK_RFCOUNT.out.rci
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

    // Enforce rf-norm pairing rules explicitly:
    // - treated may run on its own
    // - untreated requires a matching treated sample
    // - denatured requires matching treated and untreated samples
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

    ch_norm_input = ch_treated
        .join(ch_group_meta)
        .join(ch_untreated, remainder: true)
        .join(ch_denatured, remainder: true)
        .join(ch_group_rci, remainder: true)
        .map { group, treated_rcs, base_meta, untreated_rc, denatured_rc, rci_files ->
            def hasUntreated = untreated_rc ? true : false
            def hasDenatured = denatured_rc ? true : false
            def principle = (base_meta.principle ?: '').toLowerCase()
            def scoringMethod = principle == 'map'
                ? (hasUntreated ? 3 : 4)
                : (hasUntreated ? 1 : 2)
            def normMethod = scoringMethod == 2 ? 2 : 3
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

    RNAFRAMEWORK_RFFOLD (
        ch_fold_input
    )

    // Add RNAframework outputs to MultiQC input collection.
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFCOUNT.out.rc.collect { rc_file -> rc_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFCOUNT.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.xml.collect { xml_file -> xml_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFFOLD.out.structures.collect { fold_dir -> fold_dir[1] })

    // Add RNAframework versions to the MultiQC software-versions input.
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(RNAFRAMEWORK_RFCOUNT.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())

    //
    // MODULE: multiqc — aggregate pipeline quality control reports
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    summary_params      = filterSummaryParams(summary_params)
    summary_params      = addModuleOptionsSummary(summary_params, params)
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    def ch_versions_for_multiqc_yaml = softwareVersionsToYAML(ch_versions_for_multiqc_files)
        .mix(
            ch_versions_for_multiqc_tuples
                .map { process, tool, version ->
                    [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
                }
                .groupTuple(by:0)
                .map { process, tool_versions ->
                    tool_versions.unique().sort()
                    "${process}:\n${tool_versions.join('\n')}"
                }
        )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'rnastructurome_software_' + 'mqc_' + 'versions_for_multiqc.yml',
            sort: true,
            newLine: true
        )

    ch_multiqc_files = ch_multiqc_files.mix(ch_versions_for_multiqc_yaml)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    //
    // Collect software versions from all modules.
    // Modules using `path "versions.yml"` (old pattern) must be mixed in explicitly.
    // Modules using `topic: versions` (new pattern) are collected automatically by
    // channel.topic("versions") below and do NOT need to be listed here.
    //   topic-pattern: CUTADAPT_*, BOWTIE2_*, UMITOOLS_*, all SAMTOOLS_* modules
    //   file-pattern:  FASTQC, BOWTIE_BUILD/ALIGN, local modules
    //
    ch_versions = ch_versions.mix(FASTQC_PRE.out.versions.first())
    ch_versions = ch_versions.mix(FASTQC_POST.out.versions.first())
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFCOUNT.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())

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
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnastructurome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    multiqc_report   = MULTIQC.out.report.toList()        // channel: /path/to/multiqc_report.html
    mapped_bam       = ch_dedup_bam                       // channel: [ val(meta), path(bam) ]
    normalized_xml   = RNAFRAMEWORK_RFNORM.out.xml        // channel: [ val(meta), path(xml) ]
    fold_structures  = RNAFRAMEWORK_RFFOLD.out.structures // channel: [ val(meta), path(dir) ]
    versions         = ch_versions                        // channel: [ path(versions.yml) ]

}

def normaliseEnsemblSpecies(value) {
    value
        ?.toString()
        ?.trim()
        ?.toLowerCase()
        ?.replaceAll(/\s+/, '_')
}

def resolveReferenceKey(meta, fallbackOrganism) {
    def rawReference = (meta.organism ?: fallbackOrganism)?.toString()?.trim()
    if (!rawReference) {
        error("Missing organism for sample '${meta.id}'. Set organism in the samplesheet or provide --organism.")
    }
    if (rawReference.contains(' ')) {
        return normaliseEnsemblSpecies(rawReference)
    }
    if (rawReference ==~ /[a-z]+_[a-z0-9_]+/) {
        return rawReference.toLowerCase()
    }
    rawReference
}

def parseFlagstatMappedReads(flagstatFile) {
    def mappedLine = flagstatFile.readLines().find { line ->
        line ==~ /^\d+\s+\+\s+\d+\s+mapped\s+\(.*/
    }
    if (!mappedLine) {
        error("Could not parse mapped read count from flagstat file: ${flagstatFile}")
    }
    (mappedLine.tokenize()[0]) as long
}

def parseRfcountCoveredTranscripts(summaryFile) {
    def summaryLines = summaryFile.readLines().findAll { line -> line?.trim() }
    if (summaryLines.size() < 2) {
        error("Could not parse rf-count summary TSV: ${summaryFile}")
    }
    def fields = summaryLines[1].split('\t')
    if (fields.size() < 2) {
        error("rf-count summary TSV is missing the covered transcript column: ${summaryFile}")
    }
    (fields[1]) as long
}

def filterSummaryParams(summaryParams) {
    def hiddenKeys = [
        'ensembl_species_map',
        'container',
        'configFiles',
        'launchDir',
        'projectDir',
        'userName',
        'workDir'
    ] as Set

    summaryParams.collectEntries { sectionName, sectionParams ->
        if (!(sectionParams instanceof Map)) {
            return [(sectionName): sectionParams]
        }

        def filteredSection = sectionParams.findAll { key, value ->
            if (hiddenKeys.contains(key)) {
                return false
            }
            if (key == 'genomes' && value instanceof Map && value.isEmpty()) {
                return false
            }
            true
        }

        [(sectionName): filteredSection]
    }.findAll { _sectionName, sectionParams ->
        !(sectionParams instanceof Map) || !sectionParams.isEmpty()
    }
}

def addModuleOptionsSummary(summaryParams, params) {
    def sampleMetadata = parseInputSamplesheetMetadata(params.input)
    def moduleOptions = buildModuleOptionsSummary(params, sampleMetadata)
    if (moduleOptions.isEmpty()) {
        return summaryParams
    }
    summaryParams + ['Module options': moduleOptions]
}

def parseInputSamplesheetMetadata(inputPath) {
    if (!inputPath) {
        return [principles: [], conditions: [], methods: []]
    }

    def inputFile = file(inputPath.toString())
    if (!inputFile.exists()) {
        return [principles: [], conditions: [], methods: []]
    }

    def lines = inputFile.readLines().findAll { line -> line?.trim() }
    if (lines.size() < 2) {
        return [principles: [], conditions: [], methods: []]
    }

    def header = lines[0].split(',', -1)*.trim()
    def principleIdx = header.indexOf('principle')
    def conditionIdx = header.indexOf('condition')
    def methodIdx = header.indexOf('method')

    def principles = []
    def conditions = []
    def methods = []

    lines.drop(1).each { line ->
        def fields = line.split(',', -1)
        if (principleIdx >= 0 && principleIdx < fields.size()) {
            def value = fields[principleIdx]?.trim()
            if (value) principles << value
        }
        if (conditionIdx >= 0 && conditionIdx < fields.size()) {
            def value = fields[conditionIdx]?.trim()
            if (value) conditions << value
        }
        if (methodIdx >= 0 && methodIdx < fields.size()) {
            def value = fields[methodIdx]?.trim()
            if (value) methods << value
        }
    }

    [
        principles: principles.unique(),
        conditions: conditions.collect { it.toLowerCase() }.unique(),
        methods   : methods.unique()
    ]
}

def buildModuleOptionsSummary(params, sampleMetadata) {
    def moduleOptions = [:]
    def principles = (sampleMetadata.principles ?: []).collect { it.toLowerCase() }

    if (!principles || principles.contains('rt-stop')) {
        moduleOptions['bowtie_rtstop_aligner'] = 'bowtie'
        moduleOptions['bowtie_rtstop_args'] = renderBowtie1Args(params)
    }
    if (principles.contains('map')) {
        moduleOptions['bowtie_map_aligner'] = 'bowtie2'
        moduleOptions['bowtie_map_args'] = renderBowtie2Args(params)
    }

    def rfnormSummary = renderRfNormSummary(params, sampleMetadata)
    moduleOptions.putAll(rfnormSummary)

    moduleOptions.findAll { _k, v -> v != null && v.toString().trim() }
}

def renderBowtie1Args(params) {
    def manualOnly = params.bowtie_manual_only as Boolean ?: false
    def manualParams = (params.bowtie_mapping_params ?: '').toString().trim()
    if (manualOnly) {
        return manualParams ?: 'none'
    }
    def args = []
    if (params.bowtie_all as Boolean) {
        args << '-a'
    } else if (params.bowtie_k != null) {
        args << "-k ${params.bowtie_k as Integer}"
    }
    if (params.bowtie_norc as Boolean) {
        args << '--norc'
    }
    if ((params.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${params.bowtie_trim5 as Integer}"
    }
    if ((params.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${params.bowtie_trim3 as Integer}"
    }
    def seedlen = params.bowtie_seedlen != null ? params.bowtie_seedlen as Integer : 28
    args << "-l ${seedlen}"
    if (params.bowtie_v != null) {
        args << "-v ${params.bowtie_v as Integer}"
    } else {
        args << "-n ${params.bowtie_n as Integer}"
    }
    if (!(params.bowtie_all as Boolean) && params.bowtie_k == null && params.bowtie_max != null) {
        args << "-m ${params.bowtie_max as Integer}"
    }
    args << "--chunkmbs ${params.bowtie_chunkmbs as Integer}"
    if (manualParams) {
        args << manualParams
    }
    args.join(' ').trim()
}

def renderBowtie2Args(params) {
    def manualOnly = params.bowtie_manual_only as Boolean ?: false
    def manualParams = (params.bowtie_mapping_params ?: '').toString().trim()
    if (manualOnly) {
        return manualParams ?: 'none'
    }
    def args = []
    if (params.bowtie_all as Boolean) {
        args << '-a'
    } else if (params.bowtie_k != null) {
        args << "-k ${params.bowtie_k as Integer}"
    }
    if (params.bowtie_norc as Boolean) {
        args << '--norc'
    }
    if ((params.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${params.bowtie_trim5 as Integer}"
    }
    if ((params.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${params.bowtie_trim3 as Integer}"
    }
    def seedlen = params.bowtie_seedlen != null ? params.bowtie_seedlen as Integer : 22
    args << "-L ${seedlen}"
    args << "-N ${params.bowtie2_N as Integer}"
    args << "-D ${params.bowtie2_D as Integer}"
    args << "-R ${params.bowtie2_R as Integer}"
    args << "--mp ${params.bowtie2_mp}"
    args << "--dpad ${params.bowtie2_dpad as Integer}"
    args << "--rdg ${params.bowtie2_rdg}"
    args << "--rfg ${params.bowtie2_rfg}"
    if (params.bowtie2_softclip as Boolean) {
        args << '--local'
        args << "--ma ${params.bowtie2_ma as Integer}"
    }
    if (params.bowtie2_dovetail as Boolean) {
        args << '--dovetail'
    }
    if (manualParams) {
        args << manualParams
    }
    args.join(' ').trim()
}

def renderRfNormSummary(params, sampleMetadata) {
    def principles = (sampleMetadata.principles ?: []).collect { it.toLowerCase() }.unique()
    def conditions = (sampleMetadata.conditions ?: []).collect { it.toLowerCase() }.unique()
    if (principles.size() != 1) {
        return [
            rfnorm_mode: 'dynamic (mixed principles across samples)'
        ]
    }

    def principle = principles[0]
    def hasUntreated = conditions.contains('untreated')
    def hasDenatured = conditions.contains('denatured')
    def scoringMethod = principle == 'map' ? (hasUntreated ? 3 : 4) : (hasUntreated ? 1 : 2)
    def normMethod = scoringMethod == 2 ? 2 : 3
    def reactiveBases = params.rfnorm_reactive_bases ?: ((((sampleMetadata.methods ?: []).collect { it.toLowerCase() }.unique() == ['dms']) ? 'AC' : null))

    def args = [
        "-sm ${scoringMethod}",
        "-nm ${normMethod}"
    ]
    if (params.rfnorm_remap_reactivities as Boolean) args << '--remap-reactivities'
    if (reactiveBases) args << "--reactive-bases ${reactiveBases}"
    if (params.rfnorm_norm_window != null) args << "--norm-window ${params.rfnorm_norm_window as Integer}"
    if (params.rfnorm_window_offset != null) args << "--window-offset ${params.rfnorm_window_offset as Integer}"
    if (params.rfnorm_dynamic_window != null) args << "--dynamic-window ${params.rfnorm_dynamic_window as Integer}"
    if (params.rfnorm_norm_independent as Boolean) args << '--norm-independent'
    if (params.rfnorm_norm_factor) args << "--norm-factor ${params.rfnorm_norm_factor}"
    if (params.rfnorm_raw as Boolean) args << '--raw'
    if (params.rfnorm_pseudocount != null) args << "--pseudocount ${params.rfnorm_pseudocount}"
    if (params.rfnorm_max_score != null) args << "--max-score ${params.rfnorm_max_score}"
    if (params.rfnorm_ignore_lower_than_untreated as Boolean) args << '--ignore-lower-than-untreated'
    if (params.rfnorm_max_untreated_mut != null) args << "--max-untreated-mut ${params.rfnorm_max_untreated_mut}"
    if (params.rfnorm_max_mutation_rate != null) args << "--max-mutation-rate ${params.rfnorm_max_mutation_rate}"
    def meanCoverage = params.rfnorm_mean_coverage != null ? params.rfnorm_mean_coverage as BigDecimal : 0
    if (meanCoverage > 0) args << "--mean-coverage ${params.rfnorm_mean_coverage}"
    def medianCoverage = params.rfnorm_median_coverage != null ? params.rfnorm_median_coverage as BigDecimal : 0
    if (medianCoverage > 0) args << "--median-coverage ${params.rfnorm_median_coverage}"
    def nanThreshold = params.rfnorm_nan != null ? params.rfnorm_nan as Integer : 10
    if (nanThreshold != 10) args << "--nan ${nanThreshold}"
    args << '--img'
    args << "-R ${params.rnaframework_r_path}"

    [
        rfnorm_mode         : "${principle.toUpperCase()} ${hasUntreated ? 'with untreated' : 'treated-only'}${hasDenatured ? ' + denatured' : ''}",
        rfnorm_scoring      : "${rfNormScoringLabel(scoringMethod)} (sm=${scoringMethod})",
        rfnorm_normalisation: "${rfNormNormLabel(normMethod)} (nm=${normMethod})",
        rfnorm_args         : args.join(' ').trim()
    ]
}

def rfNormScoringLabel(code) {
    switch(code as Integer) {
        case 1: return 'Ding'
        case 2: return 'Rouskin'
        case 3: return 'Siegfried'
        case 4: return 'Zubradt'
        default: return "unknown"
    }
}

def rfNormNormLabel(code) {
    switch(code as Integer) {
        case 2: return '90% Winsorizing'
        case 3: return 'Box-plot normalisation'
        default: return "unknown"
    }
}

def parseCutadaptCommandArg(logFile, optionName) {
    def commandLine = logFile.readLines().find { line -> line.startsWith('Command line parameters:') }
    if (!commandLine) {
        return 'none'
    }
    def matcher = (commandLine =~ /(?:^|\s)${java.util.regex.Pattern.quote(optionName)}\s+(\S+)/)
    matcher.find() ? matcher.group(1) : 'none'
}

def countProgressionMultiqc(rows) {
    def rowEntries
    if (rows instanceof Map) {
        rowEntries = rows.entrySet().collect { [it.key, it.value] }
    } else if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else {
        rowEntries = rows
    }
    def orderedRows = rowEntries.sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sample_id = row[0]
        def metrics = row[1]
        def metricLines = metrics.collect { key, value ->
            def rendered = value instanceof BigDecimal ? String.format(java.util.Locale.ROOT, '%.2f', value) : value.toString()
            "    ${key}: ${rendered}"
        }.join('\n')
        "  ${sample_id}:\n${metricLines}"
    }.join('\n')

    """id: 'nf-core-rnastructurome-count-progression'
section_name: 'nf-core/rnastructurome Count Progression'
description: 'Mapped read retention through deduplication together with rf-count covered transcript totals.'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-count-progression'
  title: 'nf-core/rnastructurome Count Progression'
headers:
  mapped_reads_pre_dedup:
    title: 'Mapped Reads Pre-dedup'
    format: '{:,.0f}'
  mapped_reads_post_dedup:
    title: 'Mapped Reads Post-dedup'
    format: '{:,.0f}'
  reads_removed_by_dedup:
    title: 'Reads Removed by Dedup'
    format: '{:,.0f}'
  pct_removed_by_dedup:
    title: 'Dedup Removed %'
    format: '{:,.2f}'
    suffix: '%'
  rfcount_covered_transcripts:
    title: 'RFCOUNT Covered Transcripts'
    format: '{:,.0f}'
data:
${dataBlock}
"""
}

def cutadaptAdaptersMultiqc(rows) {
    def rowEntries
    if (rows instanceof Map) {
        rowEntries = rows.entrySet().collect { [it.key, it.value] }
    } else if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else {
        rowEntries = rows
    }
    def orderedRows = rowEntries.sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sample_id = row[0]
        def metrics = row[1]
        def metricLines = metrics.collect { key, value ->
            "    ${key}: '${value.toString().replace(\"'\", \"''\")}'"
        }.join('\n')
        "  ${sample_id}:\n${metricLines}"
    }.join('\n')

    """id: 'nf-core-rnastructurome-cutadapt-adapters'
section_name: 'nf-core/rnastructurome Cutadapt Adapters'
description: 'Cutadapt adapter sequences used for trimming in each sample.'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-cutadapt-adapters'
  title: 'nf-core/rnastructurome Cutadapt Adapters'
headers:
  cutadapt_mode:
    title: 'Trim Mode'
  adapter_5p:
    title: \"5' Adapter\"
  adapter_3p:
    title: \"3' Adapter\"
data:
${dataBlock}
"""
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
