//
// ALIGN_READS — align trimmed reads, coordinate-sort, deduplicate and collect
// alignment QC + library strandedness.
//
// Genome route: STAR (RT-stop SE and MaP PE). MaP PE BAMs go name-sort ->
// fixmate -> coord-sort to carry the MC tag markdup needs; RT-stop SE skips it.
// Transcriptome route: Bowtie (RT-stop) / Bowtie2 (MaP). UMI samples are
// deduplicated with umi_tools, others with samtools markdup. RSeQC infers
// per-sample strandedness from BED12 (genome route only).
//

include { STAR_ALIGN as STAR_ALIGN_RTSTOP        } from '../../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_MAP           } from '../../../modules/nf-core/star/align/main'
include { BOWTIE_ALIGN          } from '../../../modules/nf-core/bowtie/align/main'
include { BOWTIE2_ALIGN         } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT                              } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME        } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FIXMATE      } from '../../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SORT  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MARKDUP      } from '../../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_STATS        } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT     } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_PRE } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS     } from '../../../modules/nf-core/samtools/idxstats/main'
include { UMITOOLS_DEDUP        } from '../../../modules/nf-core/umitools/dedup/main'
include { BEDOPS_GTF2BED        } from '../../../modules/nf-core/bedops/gtf2bed/main'
include { RSEQC_INFEREXPERIMENT } from '../../../modules/nf-core/rseqc/inferexperiment/main'

include { buildStarAlignInputs } from '../../../workflows/rnastructurome_functions.nf'
include { resolveReferenceKey  } from '../../../workflows/rnastructurome_functions.nf'
include { collectToMap         } from '../../../workflows/rnastructurome_functions.nf'
include { parseInferExperiment } from '../../../workflows/rnastructurome_functions.nf'

workflow ALIGN_READS {

    take:
    ch_rtstop_trimmed_for_align // channel: [ val(meta), [ reads ] ]
    ch_map_trimmed_for_align    // channel: [ val(meta), [ reads ] ]
    ch_star_index               // channel: [ val(meta), path(index) ]
    ch_all_reference_gtf        // channel: [ val(meta), path(gtf) ]
    ch_bowtie_index_map         // value:   map ref_key -> [meta, index]
    ch_bowtie2_index_map        // value:   map ref_key -> [meta, index]
    ch_reference_fasta_map      // value:   map ref_key -> [meta, fasta]
    pipeline_config             // map

    main:
    ch_multiqc_files = channel.empty()

    //
    // RT-STOP ALIGNMENT
    //
    def ch_rtstop_aligned_bam      = channel.empty()

    if (!pipeline_config.transcriptome) {
        def ch_rtstop_star_split = buildStarAlignInputs(
            ch_rtstop_trimmed_for_align, ch_star_index, ch_all_reference_gtf, pipeline_config
        ).multiMap { entry ->
            reads:      entry[0]
            index:      entry[1]
            gtf:        entry[2]
            ignore_gtf: entry[3]
        }
        STAR_ALIGN_RTSTOP(
            ch_rtstop_star_split.reads,
            ch_rtstop_star_split.index,
            ch_rtstop_star_split.gtf,
            ch_rtstop_star_split.ignore_gtf
        )
        ch_rtstop_aligned_bam    = STAR_ALIGN_RTSTOP.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_RTSTOP.out.log_final.collect { _meta, log -> log })
    } else {
        def ch_rtstop_bowtie_inputs = ch_rtstop_trimmed_for_align
            .combine(ch_bowtie_index_map)
            .map { combined ->
                def meta      = combined[0]
                def reads     = combined[1]
                def index_map = combined[2]
                def ref_key   = resolveReferenceKey(meta, pipeline_config.organism)
                def idx_tuple = index_map[ref_key]
                if (!idx_tuple) error("No Bowtie index resolved for reference '${ref_key}'.")
                [ [meta, reads], idx_tuple ]
            }
        def ch_rtstop_bowtie_split = ch_rtstop_bowtie_inputs.multiMap { entry ->
            reads: entry[0]
            index: entry[1]
        }
        BOWTIE_ALIGN(ch_rtstop_bowtie_split.reads, ch_rtstop_bowtie_split.index, false)
        ch_rtstop_aligned_bam = BOWTIE_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect { _meta, log -> log })
    }

    //
    // MAP ALIGNMENT
    //
    def ch_map_aligned_bam    = channel.empty()

    if (!pipeline_config.transcriptome) {
        def ch_map_star_split = buildStarAlignInputs(
            ch_map_trimmed_for_align, ch_star_index, ch_all_reference_gtf, pipeline_config
        ).multiMap { entry ->
            reads:      entry[0]
            index:      entry[1]
            gtf:        entry[2]
            ignore_gtf: entry[3]
        }
        STAR_ALIGN_MAP(
            ch_map_star_split.reads,
            ch_map_star_split.index,
            ch_map_star_split.gtf,
            ch_map_star_split.ignore_gtf
        )
        ch_map_aligned_bam    = STAR_ALIGN_MAP.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_MAP.out.log_final.collect { _meta, log -> log })
    } else {
        def ch_map_bowtie2_inputs = ch_map_trimmed_for_align
            .combine(ch_bowtie2_index_map)
            .combine(ch_reference_fasta_map)
            .map { combined ->
                def meta      = combined[0]
                def reads     = combined[1]
                def index_map = combined[2]
                def ref_map   = combined[3]
                def ref_key   = resolveReferenceKey(meta, pipeline_config.organism)
                def idx_tuple = index_map[ref_key]
                def fasta_t   = ref_map[ref_key]
                if (!idx_tuple) error("No Bowtie2 index resolved for reference '${ref_key}'.")
                if (!fasta_t)   error("No transcript FASTA resolved for reference '${ref_key}' in MaP alignment.")
                [ [meta, reads], idx_tuple, [ fasta_t[0], fasta_t[1] ] ]
            }
        def ch_map_bowtie2_split = ch_map_bowtie2_inputs.multiMap { entry ->
            reads: entry[0]
            index: entry[1]
            fasta: entry[2]
        }
        BOWTIE2_ALIGN(ch_map_bowtie2_split.reads, ch_map_bowtie2_split.index, ch_map_bowtie2_split.fasta, false, false)
        ch_map_aligned_bam = BOWTIE2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { _meta, log -> log })
    }

    //
    // MaP (PE) BAMs need name-sort → fixmate → coordinate-sort to add the MC tag
    // required by samtools markdup. RT-stop (SE) BAMs skip this step.
    //
    SAMTOOLS_SORT_NAME (
        ch_map_aligned_bam,
        channel.value([ [], [], [] ]),
        false
    )
    SAMTOOLS_FIXMATE (
        SAMTOOLS_SORT_NAME.out.bam
    )
    ch_mapped_bam = ch_rtstop_aligned_bam.mix(SAMTOOLS_FIXMATE.out.bam)

    //
    // MODULE: samtools sort — coordinate-sort mapped (genome) BAMs
    // FASTA/FAI not needed for BAM output (only required for CRAM); pass empty.
    //
    SAMTOOLS_SORT (
        ch_mapped_bam,
        channel.value([ [], [], [] ]),
        false
    )

    //
    // MODULE: samtools index (sorted) — index sorted genome BAMs
    //
    SAMTOOLS_INDEX_SORT (
        SAMTOOLS_SORT.out.bam
    )

    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_SORT.out.index.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
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
    // MODULE: samtools markdup — deduplicate non-UMI BAMs
    // FASTA/FAI not needed for BAM output; pass empty.
    // Drop the bai from [meta, bam, bai] — markdup only takes [meta, bam].
    //
    def ch_non_umi_bam = dedup_branches.non_umi.map { meta, bam, _bai -> [ meta, bam ] }
    if (!params.skip_markdup) {
        SAMTOOLS_MARKDUP (
            ch_non_umi_bam,
            channel.value([ [], [], [] ])
        )
        ch_dedup_bam = UMITOOLS_DEDUP.out.bam.mix(SAMTOOLS_MARKDUP.out.bam)
    } else {
        ch_dedup_bam = UMITOOLS_DEDUP.out.bam.mix(ch_non_umi_bam)
    }

    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_DEDUP.out.log.collect { dedup_log -> dedup_log[1] })

    //
    // MODULE: samtools index (final) — index deduplicated BAMs
    //
    SAMTOOLS_INDEX_FINAL (
        ch_dedup_bam
    )

    ch_markdup_bam_bai = ch_dedup_bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_FINAL.out.index.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
        .map { _sample_id, bam_tuple, bai_tuple ->
            [ bam_tuple[0], bam_tuple[1], bai_tuple[1] ]
        }

    //
    // MODULE: samtools stats — collect alignment statistics
    // FASTA/FAI not needed for transcript BAMs; pass empty.
    //
    SAMTOOLS_STATS (
        ch_markdup_bam_bai,
        channel.value([ [], [], [] ])
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

    // Software versions are collected later via ch_versions (old-style emit: versions)
    // and channel.topic("versions") (new-style topic-based modules)

    //
    // MODULES: BEDOPS_GTF2BED + RSEQC_INFEREXPERIMENT (STAR route only)
    // Convert each reference GTF to BED12 once, then run infer_experiment on
    // every final BAM to determine library strandedness automatically.
    //
    def ch_strandedness_by_id = channel.empty()

    if (!pipeline_config.transcriptome) {
        BEDOPS_GTF2BED(ch_all_reference_gtf)

        def ch_reference_bed_map = collectToMap(
            BEDOPS_GTF2BED.out.bed.map { meta, bed -> [ meta.id.toString(), [meta, bed] ] }
        )

        def ch_infer_inputs = ch_markdup_bam_bai
            .combine(ch_reference_bed_map)
            .map { meta, bam, bai, bed_map ->
                def ref_key = resolveReferenceKey(meta, pipeline_config.organism)
                def bed_t   = bed_map[ref_key]
                // Skip samples whose reference has no usable BED (e.g. viral synthetic GTFs)
                bed_t ? [ [meta, bam, bai], bed_t[1] ] : null
            }
            .filter { entry -> entry != null }

        def ch_infer_split = ch_infer_inputs.multiMap { entry ->
            bam: entry[0]
            bed: entry[1]
        }

        RSEQC_INFEREXPERIMENT(ch_infer_split.bam, ch_infer_split.bed)
        ch_multiqc_files = ch_multiqc_files.mix(
            RSEQC_INFEREXPERIMENT.out.txt.collect { _meta, txt -> txt }
        )

        ch_strandedness_by_id = RSEQC_INFEREXPERIMENT.out.txt
            .map { meta, txt -> [ meta.id.toString(), parseInferExperiment(txt) ] }
    }

    emit:
    markdup_bam_bai    = ch_markdup_bam_bai               // channel: [ val(meta), path(bam), path(bai) ]
    dedup_bam          = ch_dedup_bam                     // channel: [ val(meta), path(bam) ]
    strandedness_by_id = ch_strandedness_by_id            // channel: [ id, strandedness ]
    flagstat_pre       = SAMTOOLS_FLAGSTAT_PRE.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    flagstat_post      = SAMTOOLS_FLAGSTAT.out.flagstat      // channel: [ val(meta), path(flagstat) ]
    multiqc_files      = ch_multiqc_files
}
