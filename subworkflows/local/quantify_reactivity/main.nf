//
// QUANTIFY_REACTIVITY — per-base RT-stop / mutation counts as RNAFramework RC files.
//
// Genome route (STAR): rf-count-genome on the genome BAM (multi-mappers counted
// at all genome positions), then rf-rctools extract redistributes counts to
// transcripts via the reference GTF — more accurate than STAR TranscriptomeSAM.
// Transcriptome route (--transcriptome): rf-count directly on the transcript BAM.
//

include { RNAFRAMEWORK_RFCOUNT           } from '../../../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFCOUNT_GENOME    } from '../../../modules/local/rnaframework/count_genome/main'
include { RNAFRAMEWORK_RFRCTOOLS_EXTRACT } from '../../../modules/local/rnaframework/rctools/extract/main'

include { collectToMap        } from '../../../workflows/rnastructurome_functions.nf'
include { resolveReferenceKey } from '../../../workflows/rnastructurome_functions.nf'

workflow QUANTIFY_REACTIVITY {

    take:
    ch_markdup_bam_bai             // channel: [ val(meta), path(bam), path(bai) ]
    ch_strandedness_by_id          // channel: [ id, strandedness ]
    ch_reference_genome_fasta_keyed // channel: [ key, [meta, fasta] ]
    ch_reference_gtf_map           // value:   map ref_key -> [meta, gtf]
    ch_reference_fasta_map         // value:   map ref_key -> [meta, fasta]
    pipeline_config                // map

    main:
    ch_versions = channel.empty()

    // MODULE: rf-count / rf-count-genome — per-base RT-stop or mutation counts
    // Genome route (STAR): rf-count-genome with genome FASTA + genome-coord BAM
    // Transcriptome route (--transcriptome): rf-count with transcript FASTA
    //
    def ch_rfcount_rc        = channel.empty()
    def ch_rfcount_rci       = channel.empty()
    def ch_rfcount_summary   = channel.empty()
    def ch_rfcount_plots     = channel.empty()

    //
    // MODULES: rf-count — per-base RT-stop or mutation counts
    //
    // STAR route (genome alignment): genome BAM → rf-count-genome → rf-rctools extract → transcript RC files
    //   Multi-mappers are counted at all genome positions then redistributed to transcripts via GTF,
    //   which is more accurate than STAR's internal TranscriptomeSAM quantification.
    //
    // Bowtie route (--transcriptome): transcript-coordinate BAM → rf-count → transcript RC files
    //
    if (!pipeline_config.transcriptome) {
        def ch_genome_fasta_map = collectToMap(ch_reference_genome_fasta_keyed)

        // Annotate each BAM with the per-sample strandedness inferred by RSeQC.
        // remainder: true keeps samples whose reference had no BED (e.g. viral) — they
        // get strandedness = null, and modules.config falls back to params.rfcount_strandedness.
        def ch_bam_stranded = ch_markdup_bam_bai
            .map { meta, bam, bai -> [ meta.id.toString(), meta, bam, bai ] }
            .join(ch_strandedness_by_id, remainder: true)
            .map { _id, meta, bam, bai, strandedness ->
                [ meta + [strandedness: strandedness], bam, bai ]
            }

        def ch_rfcount_genome_inputs = ch_bam_stranded
            .combine(ch_genome_fasta_map)
            .map { combined ->
                def meta    = combined[0]
                def bam     = combined[1]
                def bai     = combined[2]
                def ref_map = combined[3]
                def ref_key = resolveReferenceKey(meta, pipeline_config.organism)
                def fasta_t = ref_map[ref_key]
                if (!fasta_t) error("No genome FASTA resolved for reference '${ref_key}' for rf-count-genome.")
                [ [meta, bam, bai], fasta_t ]
            }
        def ch_rg_split = ch_rfcount_genome_inputs.multiMap { entry ->
            bam:   entry[0]
            fasta: entry[1]
        }
        RNAFRAMEWORK_RFCOUNT_GENOME(ch_rg_split.bam, ch_rg_split.fasta)
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT_GENOME.out.plots
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFCOUNT_GENOME.out.versions)

        // rf-rctools extract: genome RC → transcript-level RC using reference GTF.
        // The module generates per-file .rci indexes itself (rf-rctools index) and calls
        // rf-rctools extract with the BASENAME so strand-aware extraction is activated.
        // It also rewrites the rf-count-genome summary's 'covered' (a genome reference count)
        // with the covered-transcript count, so pair the RC with its summary (same task, strict
        // join — both always emit) and pass it through.
        def ch_rc_with_summary = RNAFRAMEWORK_RFCOUNT_GENOME.out.rc
            .map { meta, rc -> [ meta.id.toString(), meta, rc ] }
            .join(RNAFRAMEWORK_RFCOUNT_GENOME.out.summary.map { meta, summary -> [ meta.id.toString(), summary ] })
            .map { _id, meta, rc, summary -> [ meta, rc, summary ] }

        def ch_rctools_inputs = ch_rc_with_summary
            .combine(ch_reference_gtf_map)
            .map { combined ->
                def meta    = combined[0]
                def rc      = combined[1]
                def summary = combined[2]
                def gtf_map = combined[3]
                def ref_key = resolveReferenceKey(meta, pipeline_config.organism)
                def gtf_t   = gtf_map[ref_key]
                if (!gtf_t) error("No GTF resolved for reference '${ref_key}' for rf-rctools extract.")
                [ [meta, rc, [], summary], gtf_t ]
            }
        def ch_rct_split = ch_rctools_inputs.multiMap { entry ->
            rc:  entry[0]
            gtf: entry[1]
        }
        RNAFRAMEWORK_RFRCTOOLS_EXTRACT(ch_rct_split.rc, ch_rct_split.gtf)
        ch_rfcount_rc      = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rc
        ch_rfcount_rci     = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rci
        ch_rfcount_summary = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.summary
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.versions)
    } else {
        def ch_rfcount_inputs = ch_markdup_bam_bai
            .combine(ch_reference_fasta_map)
            .map { combined ->
                def meta    = combined[0]
                def bam     = combined[1]
                def bai     = combined[2]
                def ref_map = combined[3]
                def ref_key = resolveReferenceKey(meta, pipeline_config.organism)
                def fasta_t = ref_map[ref_key]
                if (!fasta_t) error("No transcript FASTA resolved for reference '${ref_key}' for rf-count.")
                [ [meta, bam, bai], fasta_t ]
            }
        def ch_rc_split = ch_rfcount_inputs.multiMap { entry ->
            bam:   entry[0]
            fasta: entry[1]
        }
        RNAFRAMEWORK_RFCOUNT(ch_rc_split.bam, ch_rc_split.fasta)
        ch_rfcount_rc      = RNAFRAMEWORK_RFCOUNT.out.rc
        ch_rfcount_rci     = RNAFRAMEWORK_RFCOUNT.out.rci
        ch_rfcount_summary = RNAFRAMEWORK_RFCOUNT.out.summary
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT.out.plots
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFCOUNT.out.versions)
    }

    emit:
    rc       = ch_rfcount_rc       // channel: [ val(meta), path(rc) ]
    rci      = ch_rfcount_rci      // channel: [ val(meta), path(rci) ]
    summary  = ch_rfcount_summary  // channel: [ val(meta), path(summary_tsv) ]
    plots    = ch_rfcount_plots    // channel: [ val(meta), path(plots) ]
    versions = ch_versions
}
