// Subworkflow with functionality specific to the nf-core/rnastructurome pipeline

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    // Print version and exit if required and dump pipeline parameters to JSON file
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    // Validate parameters and generate parameter summary to stdout
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/rnastructurome ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/rnastructurome/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    // Check config provided to the pipeline
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    // Create channel from input file provided through `input`

    def samplesheet_rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    warnOnSamplesheetOverrides(samplesheet_rows)

    channel
        .fromList(samplesheet_rows)
        .map {
            meta, fastq_1, fastq_2 ->
                def resolved_meta = meta + [
                    sample_id     : meta.sample_id ?: params.sample_id,
                    library_layout: fastq_2 ? 'PAIRED' : 'SINGLE',
                    method        : meta.method ?: params.method,
                    principle     : meta.principle ?: params.principle,
                    chemical      : meta.chemical ?: params.chemical,
                    RT_enzyme     : meta.RT_enzyme ?: params.RT_enzyme,
                    pH            : hasMetadataValue(meta.pH) ? meta.pH : params.pH,
                    organism      : meta.organism ?: params.organism,
                    adapter_3p    : meta.adapter_3p,
                    adapter_5p    : meta.adapter_5p,
                    umi_pattern   : meta.umi_pattern ?: params.umi_pattern,
                    sample_group  : meta.sample_group,
                    replicate     : meta.replicate
                ]

                if (!fastq_2) {
                    return [ resolved_meta.id, resolved_meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ resolved_meta.id, resolved_meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Completion email and summary
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def normaliseMetadataValue(value) {
    return value.toString().trim().toLowerCase()
}

// The per-sample params are fallbacks for rows leaving the column blank — a row that sets it wins.
// Warn once per field so the precedence is visible rather than silent. Compared case-insensitively,
// as principle/method are read downstream: warning about a value the user did not actually
// contradict is worse than staying quiet.
def warnOnSamplesheetOverrides(rows) {
    ['sample_id', 'method', 'principle', 'chemical', 'RT_enzyme', 'pH', 'organism', 'umi_pattern'].each { field ->
        def param_value = params[field]
        if (!hasMetadataValue(param_value)) {
            return
        }
        def conflicting = rows
            .collect { row -> row[0] }
            .findAll { meta ->
                hasMetadataValue(meta[field]) &&
                    normaliseMetadataValue(meta[field]) != normaliseMetadataValue(param_value)
            }
        if (conflicting) {
            log.warn(
                "--${field}=${param_value} is overridden by the samplesheet for ${conflicting.size()} sample(s) " +
                "(e.g. '${conflicting[0].id}' = ${conflicting[0][field]}); the samplesheet takes precedence."
            )
        }
    }
}

def hasMetadataValue(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Collection) {
        return !value.isEmpty()
    }
    return value.toString().trim()
}

// Validate channels from input samplesheet
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
// Tools credited in the MultiQC methods section, in pipeline order. Citation and reference live in one
// entry so the two lists cannot drift; `when` mirrors the conditions the workflow actually branches on.
def pipelineToolReferences(transcriptome) {
    def counts_on_genome = !transcriptome && params.count_genome
    [
        [ when: true,
          cite: 'FastQC (Andrews 2010)',
          ref : '<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>' ],
        [ when: params.umi_pattern as Boolean,
          cite: 'UMI-tools (Smith et al. 2017)',
          ref : '<li>Smith T, et al. (2017). UMI-tools: modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491–499. doi: 10.1101/gr.209601.116</li>' ],
        [ when: true,
          cite: 'Cutadapt (Martin 2011)',
          ref : '<li>Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. doi: 10.14806/ej.17.1.200</li>' ],
        [ when: !transcriptome && !params.count_genome,
          cite: 'GffRead (Pertea & Pertea 2020)',
          ref : '<li>Pertea G & Pertea M (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. doi: 10.12688/f1000research.23297.2</li>' ],
        [ when: !transcriptome,
          cite: 'STAR (Dobin et al. 2013)',
          ref : '<li>Dobin A, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. doi: 10.1093/bioinformatics/bts635</li>' ],
        [ when: transcriptome as Boolean,
          cite: 'Bowtie (Langmead et al. 2009)',
          ref : '<li>Langmead B, et al. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 10(3), R25. doi: 10.1186/gb-2009-10-3-r25</li>' ],
        [ when: transcriptome as Boolean,
          cite: 'Bowtie2 (Langmead & Salzberg 2012)',
          ref : '<li>Langmead B & Salzberg SL (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. doi: 10.1038/nmeth.1923</li>' ],
        [ when: true,
          cite: 'SAMtools (Danecek et al. 2021)',
          ref : '<li>Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>' ],
        [ when: counts_on_genome,
          cite: 'BEDOPS (Neph et al. 2012)',
          ref : '<li>Neph S, et al. (2012). BEDOPS: high-performance genomic feature operations. Bioinformatics, 28(14), 1919–1920. doi: 10.1093/bioinformatics/bts277</li>' ],
        [ when: counts_on_genome,
          cite: 'RSeQC (Wang et al. 2012)',
          ref : '<li>Wang L, Wang S & Li W (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics, 28(16), 2184–2185. doi: 10.1093/bioinformatics/bts356</li>' ],
        [ when: true,
          cite: 'RNAFramework (Incarnato et al. 2018)',
          ref : '<li>Incarnato D, et al. (2018). RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications. Nucleic Acids Research, 46(W1), W121–W127. doi: 10.1093/nar/gky486</li>' ],
        [ when: !params.stop_after_jackknife,
          cite: 'ViennaRNA (Lorenz et al. 2011)',
          ref : '<li>Lorenz R, et al. (2011). ViennaRNA Package 2.0. Algorithms for Molecular Biology, 6, 26. doi: 10.1186/1748-7188-6-26</li>' ],
        [ when: !params.stop_after_jackknife && params.r2dt,
          cite: 'R2DT (Sweeney et al. 2021)',
          ref : '<li>Sweeney BA, et al. (2021). R2DT is a framework for predicting and visualising RNA secondary structure using templates. Nature Communications, 12(1), 3494. doi: 10.1038/s41467-021-23555-5</li>' ],
        [ when: !params.stop_after_jackknife,
          cite: 'UCSC wigToBigWig (Kent et al. 2010)',
          ref : '<li>Kent WJ, et al. (2010). BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204–2207. doi: 10.1093/bioinformatics/btq351</li>' ],
        [ when: true,
          cite: 'MultiQC (Ewels et al. 2016)',
          ref : '<li>Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: 10.1093/bioinformatics/btw354</li>' ]
    ].findAll { tool -> tool.when }
}

def toolCitationText(transcriptome) {
    def tools = pipelineToolReferences(transcriptome).collect { tool -> tool.cite }
    return "Tools used in the workflow included: ${tools.join(', ')}."
}

def toolBibliographyText(transcriptome) {
    return pipelineToolReferences(transcriptome).collect { tool -> tool.ref }.join(' ')
}

def methodsDescriptionText(mqc_methods_yaml, transcriptome) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Loop to handle multiple DOIs, stripping `https://doi.org/` (DOIs vs resolvers) and spaces
        // (manifest.doi is a string, not a proper list).
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText(transcriptome)
    meta["tool_bibliography"] = toolBibliographyText(transcriptome)


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
