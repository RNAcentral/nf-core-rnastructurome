/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CUTADAPT as CUTADAPT_RTSTOP } from '../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_MAP    } from '../modules/nf-core/cutadapt/main'
include { BOWTIE_BUILD          } from '../modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN          } from '../modules/nf-core/bowtie/align/main'
include { BOWTIE2_BUILD         } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN         } from '../modules/nf-core/bowtie2/align/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
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

    // Branch by probing principle so RT-stop and MaP can use different default cutadapt args
    def principle_branches = ch_samplesheet_checked.branch { meta, reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    //
    // MODULE: Run cutadapt
    //
    CUTADAPT_RTSTOP (
        principle_branches.rtstop
    )

    CUTADAPT_MAP (
        principle_branches.map
    )

    ch_trimmed_reads = CUTADAPT_RTSTOP.out.reads.mix(CUTADAPT_MAP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_RTSTOP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_MAP.out.log.collect { cutadapt_log -> cutadapt_log[1] })

    ch_rtstop_genome_build_fasta = principle_branches.rtstop
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}'. Set genome_build in samplesheet or provide --genome_build/--genome.")
            }
            def genome_fasta = params.genomes?.containsKey(genome_build) ? params.genomes[genome_build]?.fasta : null
            def fasta_path = genome_fasta ?: params.fasta
            if (!fasta_path) {
                error("No FASTA configured for genome_build '${genome_build}'. Add params.genomes['${genome_build}'].fasta or provide --fasta.")
            }
            [ genome_build, fasta_path.toString() ]
        }
        .groupTuple()
        .map { genome_build, fasta_paths ->
            def unique_fasta_paths = fasta_paths.unique()
            if (unique_fasta_paths.size() != 1) {
                error("Multiple FASTA paths were resolved for genome_build '${genome_build}': ${unique_fasta_paths.join(', ')}")
            }
            [ [ id: genome_build, genome_build: genome_build ], file(unique_fasta_paths[0], checkIfExists: true) ]
        }

    ch_map_genome_build_fasta = principle_branches.map
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}'. Set genome_build in samplesheet or provide --genome_build/--genome.")
            }
            def genome_fasta = params.genomes?.containsKey(genome_build) ? params.genomes[genome_build]?.fasta : null
            def fasta_path = genome_fasta ?: params.fasta
            if (!fasta_path) {
                error("No FASTA configured for genome_build '${genome_build}'. Add params.genomes['${genome_build}'].fasta or provide --fasta.")
            }
            [ genome_build, fasta_path.toString() ]
        }
        .groupTuple()
        .map { genome_build, fasta_paths ->
            def unique_fasta_paths = fasta_paths.unique()
            if (unique_fasta_paths.size() != 1) {
                error("Multiple FASTA paths were resolved for genome_build '${genome_build}': ${unique_fasta_paths.join(', ')}")
            }
            [ [ id: genome_build, genome_build: genome_build ], file(unique_fasta_paths[0], checkIfExists: true) ]
        }

    //
    // MODULE: Build Bowtie/Bowtie2 indices
    //
    BOWTIE_BUILD (
        ch_rtstop_genome_build_fasta
    )

    BOWTIE2_BUILD (
        ch_map_genome_build_fasta
    )

    ch_bowtie_index_keyed = BOWTIE_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
    ch_bowtie2_index_keyed = BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
    ch_reference_fasta_keyed = ch_map_genome_build_fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }

    ch_rtstop_align_inputs = CUTADAPT_RTSTOP.out.reads
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, reads] ]
        }
        .join(ch_bowtie_index_keyed)
        .map { genome_build, reads_tuple, index_tuple -> [ reads_tuple, index_tuple ] }

    ch_map_align_inputs = CUTADAPT_MAP.out.reads
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, reads] ]
        }
        .join(ch_bowtie2_index_keyed)
        .join(ch_reference_fasta_keyed)
        .map { genome_build, reads_tuple, index_tuple, fasta_tuple -> [ reads_tuple, index_tuple, fasta_tuple ] }

    //
    // MODULE: Align RT-stop with Bowtie v1 and MaP with Bowtie2
    //
    BOWTIE_ALIGN (
        ch_rtstop_align_inputs.map { reads_tuple, index_tuple -> reads_tuple },
        ch_rtstop_align_inputs.map { reads_tuple, index_tuple -> index_tuple },
        false
    )

    BOWTIE2_ALIGN (
        ch_map_align_inputs.map { reads_tuple, index_tuple, fasta_tuple -> reads_tuple },
        ch_map_align_inputs.map { reads_tuple, index_tuple, fasta_tuple -> index_tuple },
        ch_map_align_inputs.map { reads_tuple, index_tuple, fasta_tuple -> fasta_tuple },
        false,
        false
    )

    ch_mapped_bam = BOWTIE_ALIGN.out.bam.mix(BOWTIE2_ALIGN.out.bam)
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect { bowtie_log -> bowtie_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { bowtie2_log -> bowtie2_log[1] })

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { fastqc_zip -> fastqc_zip[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

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
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
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
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
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

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    mapped_bam     = ch_mapped_bam               // channel: [ val(meta), path(bam) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
