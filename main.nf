#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnastructurome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnastructurome
    Website: https://nf-co.re/rnastructurome
    Slack  : https://nfcore.slack.com/channels/rnastructurome
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASTRUCTUROME  } from './workflows/rnastructurome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnastructurome_pipeline'
include { allReferencesUseNcbiRoute } from './workflows/rnastructurome_functions.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main analysis pipeline depending on type of input
workflow NFCORE_RNASTRUCTUROME {

    take:
    samplesheet      // channel: samplesheet read in from --input
    pipeline_config  // map: pipeline configuration captured at the entry workflow

    main:

    // WORKFLOW: Run pipeline
    RNASTRUCTUROME (
        samplesheet,
        pipeline_config
    )
    emit:
    multiqc_report  = RNASTRUCTUROME.out.multiqc_report  // channel: /path/to/multiqc_report.html
    normalized_xml  = RNASTRUCTUROME.out.normalized_xml  // channel: [ val(meta), path(xml) ]
    fold_structures = RNASTRUCTUROME.out.fold_structures // channel: [ val(meta), path(dir) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    def pipeline_config = buildPipelineConfig(params)

    // NCBI (bacteria/viral) references have no introns, so STAR offers nothing over Bowtie — auto-enable
    // the transcriptome route when every reference resolves to NCBI, so users don't need --transcriptome.
    if (!pipeline_config.transcriptome
            && allReferencesUseNcbiRoute(params.input, "${projectDir}/assets/schema_input.json", pipeline_config)) {
        log.info('[rnastructurome] All references resolve to the NCBI route (no-intron organisms) — enabling the transcriptome (Bowtie) route automatically. Pass --transcriptome to set it explicitly.')
        pipeline_config = pipeline_config + [ transcriptome: true ]
    }

    // A user-supplied --fasta without --gtf can only be used on the transcriptome route: the genome route
    // needs the GTF to extract and count transcripts. Assume the FASTA is a transcriptome and enable the
    // route automatically rather than failing deep in STAR — pass --gtf to use the genome route instead.
    if (!pipeline_config.transcriptome && pipeline_config.fasta && !pipeline_config.gtf) {
        log.warn('[rnastructurome] --fasta was supplied without --gtf; enabling the transcriptome route.')
        pipeline_config = pipeline_config + [ transcriptome: true ]
    }
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        pipeline_config
    )

    // WORKFLOW: Run main workflow
    NFCORE_RNASTRUCTUROME (
        PIPELINE_INITIALISATION.out.samplesheet,
        pipeline_config
    )
    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RNASTRUCTUROME.out.multiqc_report,
        params.max_multiqc_email_size
    )
}

def buildPipelineConfig(all_params) {
    [
        sample_id                        : all_params.sample_id,
        method                           : all_params.method,
        principle                        : all_params.principle,
        chemical                         : all_params.chemical,
        RT_enzyme                        : all_params.RT_enzyme,
        pH                               : all_params.pH,
        organism                         : all_params.organism,
        fasta                            : all_params.fasta,
        transcriptome                    : all_params.transcriptome,
        count_genome                     : all_params.count_genome,
        gtf                              : all_params.gtf,
        genomes                          : all_params.genomes,
        umi_pattern                      : all_params.umi_pattern,
        ensembl_base_url                 : all_params.ensembl_base_url,
        ensembl_release                  : all_params.ensembl_release,
        ensembl_species_map              : all_params.ensembl_species_map,
        ncbi_accessions_map              : all_params.ncbi_accessions_map,
        multiqc_config                   : all_params.multiqc_config,
        multiqc_logo                     : all_params.multiqc_logo,
        multiqc_methods_description      : all_params.multiqc_methods_description,
        outdir                           : all_params.outdir,
        input                            : all_params.input,
        bowtie_k                         : all_params.bowtie_k,
        bowtie_all                       : all_params.bowtie_all,
        bowtie_trim5                     : all_params.bowtie_trim5,
        bowtie_trim3                     : all_params.bowtie_trim3,
        bowtie_n                         : all_params.bowtie_n,
        bowtie_v                         : all_params.bowtie_v,
        bowtie_max                       : all_params.bowtie_max,
        bowtie_chunkmbs                  : all_params.bowtie_chunkmbs,
        bowtie2_preset                   : all_params.bowtie2_preset,
        bowtie2_mp                       : all_params.bowtie2_mp,
        bowtie2_dpad                     : all_params.bowtie2_dpad,
        bowtie2_rdg                      : all_params.bowtie2_rdg,
        bowtie2_rfg                      : all_params.bowtie2_rfg,
        bowtie2_softclip                 : all_params.bowtie2_softclip,
        bowtie2_ma                       : all_params.bowtie2_ma,
        bowtie2_dovetail                 : all_params.bowtie2_dovetail,
        jackknife_reference              : all_params.jackknife_reference,
        rfeval_reference                 : all_params.rfeval_reference,
        rfjackknife_pool_all             : all_params.rfjackknife_pool_all,
        stop_after_jackknife             : all_params.stop_after_jackknife,
        structextract                    : all_params.structextract,
        correlate_replicates             : all_params.correlate_replicates,
        rfnorm_use_normfactor            : all_params.rfnorm_use_normfactor,
        rfnorm_reactive_bases            : all_params.rfnorm_reactive_bases,
        rfnorm_remap_reactivities        : all_params.rfnorm_remap_reactivities,
        rfnorm_norm_window               : all_params.rfnorm_norm_window,
        rfnorm_window_offset             : all_params.rfnorm_window_offset,
        rfnorm_dynamic_window            : all_params.rfnorm_dynamic_window,
        rfnorm_norm_independent          : all_params.rfnorm_norm_independent,
        rfnorm_score_method              : all_params.rfnorm_score_method,
        rfnorm_norm_method               : all_params.rfnorm_norm_method,
        rfnorm_raw                       : all_params.rfnorm_raw,
        rfnorm_pseudocount               : all_params.rfnorm_pseudocount,
        rfnorm_ignore_lower_than_untreated: all_params.rfnorm_ignore_lower_than_untreated,
        rfnorm_mean_coverage             : all_params.rfnorm_mean_coverage,
        rfnorm_median_coverage           : all_params.rfnorm_median_coverage,
        rfnorm_nan                       : all_params.rfnorm_nan,
        rnaframework_r_path              : all_params.rnaframework_r_path
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
