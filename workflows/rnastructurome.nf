/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_PRE   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST  } from '../modules/nf-core/fastqc/main'
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
include { SAMTOOLS_IDXSTATS     } from '../modules/nf-core/samtools/idxstats/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_FAIDX        } from '../modules/nf-core/samtools/faidx/main'
include { RNAFRAMEWORK_RFCOUNT  } from '../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFNORM   } from '../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFFOLD   } from '../modules/local/rnaframework/fold/main'
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

    def ch_pretrim_fastqc_input = ch_samplesheet_checked
    def ch_samplesheet_for_branching = ch_samplesheet_checked

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

    //
    // MODULE: fastqc (post-trim) — quality control on trimmed reads
    //
    FASTQC_POST (
        ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    ch_rtstop_genome_build_fasta = principle_branches.rtstop
        .map { meta, _reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}'. Set genome_build in samplesheet or provide --genome_build/--genome.")
            }
            def genome_fasta = params.genomes?.containsKey(genome_build) ? params.genomes[genome_build]?.fasta : null
            def fasta_path = params.fasta ?: genome_fasta
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
        .map { meta, _reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}'. Set genome_build in samplesheet or provide --genome_build/--genome.")
            }
            def genome_fasta = params.genomes?.containsKey(genome_build) ? params.genomes[genome_build]?.fasta : null
            def fasta_path = params.fasta ?: genome_fasta
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

    ch_all_genome_build_fasta = ch_rtstop_genome_build_fasta
        .mix(ch_map_genome_build_fasta)
        .map { meta, fasta ->
            [ meta.id.toString(), fasta.toString() ]
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
    // MODULE: bowtie-build — build Bowtie v1 indices for RT-stop alignment
    //
    BOWTIE_BUILD (
        ch_rtstop_genome_build_fasta
    )

    //
    // MODULE: bowtie2-build — build Bowtie2 indices for MaP alignment
    //
    BOWTIE2_BUILD (
        ch_map_genome_build_fasta
    )

    ch_bowtie_index_keyed = BOWTIE_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
    ch_bowtie2_index_keyed = BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
    ch_reference_fasta_keyed = ch_all_genome_build_fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }

    ch_rtstop_align_inputs = CUTADAPT_RTSTOP.out.reads
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, reads] ]
        }
        .join(ch_bowtie_index_keyed)
        .map { _genome_build, reads_tuple, index_tuple -> [ reads_tuple, index_tuple ] }

    ch_map_align_inputs = CUTADAPT_MAP.out.reads
        .map { meta, reads ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, reads] ]
        }
        .join(ch_bowtie2_index_keyed)
        .join(ch_reference_fasta_keyed)
        .map { _genome_build, reads_tuple, index_tuple, fasta_tuple -> [ reads_tuple, index_tuple, fasta_tuple ] }

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

    //
    // MODULE: samtools faidx — index reference FASTA for markdup
    //
    SAMTOOLS_FAIDX (
        ch_all_genome_build_fasta.map { meta, fasta -> [meta, fasta, []] },
        false
    )

    ch_reference_fasta_fai_keyed = SAMTOOLS_FAIDX.out.fai
        .map { meta, fai -> [ meta.id.toString(), [meta, fai] ] }
        .join(ch_reference_fasta_keyed)
        .map { genome_build, fai_tuple, fasta_tuple ->
            [ genome_build, [fasta_tuple[0], fasta_tuple[1], fai_tuple[1]] ]
        }

    ch_sorted_inputs = ch_mapped_bam
        .map { meta, bam ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}' while preparing samtools sort input.")
            }
            [ genome_build, [meta, bam] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _genome_build, bam_tuple, fasta_tuple ->
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
    // MODULE: umi_tools dedup — deduplicate UMI-tagged BAMs
    //
    UMITOOLS_DEDUP (
        dedup_branches.umi,
        false
    )

    ch_markdup_inputs = dedup_branches.non_umi
        .map { meta, bam, _bai ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}' while preparing samtools markdup input.")
            }
            [ genome_build, [meta, bam] ]
        }
        .join(ch_reference_fasta_fai_keyed)
        .map { _genome_build, bam_tuple, fasta_fai_tuple ->
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
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}' while preparing samtools stats input.")
            }
            [ genome_build, [meta, bam, bai] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _genome_build, bam_bai_tuple, fasta_tuple ->
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
    // Collate pre-RNAFramework software versions so MultiQC can run in parallel
    // with rf-count/rf-norm/rf-fold instead of waiting for the full workflow.
    //
    ch_versions_for_multiqc_files = channel.empty()
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(FASTQC_PRE.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(FASTQC_POST.out.versions.first())
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(BOWTIE_BUILD.out.versions)
    ch_versions_for_multiqc_files = ch_versions_for_multiqc_files.mix(BOWTIE_ALIGN.out.versions)

    ch_versions_for_multiqc_tuples = channel.empty()
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(UMITOOLS_EXTRACT.out.versions_umitools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(UMITOOLS_DEDUP.out.versions_umitools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(CUTADAPT_RTSTOP.out.versions_cutadapt)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(CUTADAPT_MAP.out.versions_cutadapt)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_BUILD.out.versions_bowtie2)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_bowtie2)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(BOWTIE2_ALIGN.out.versions_pigz)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_SORT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_INDEX_SORT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_MARKDUP.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_INDEX_FINAL.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_STATS.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_FLAGSTAT.out.versions_samtools)
    ch_versions_for_multiqc_tuples = ch_versions_for_multiqc_tuples.mix(SAMTOOLS_IDXSTATS.out.versions_samtools)

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
            name: 'nf_core_' + 'rnastructurome_software_' + 'mqc_' + 'versions_pre_rnaframework.yml',
            sort: true,
            newLine: true
        )

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
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

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

    def ch_multiqc_done = MULTIQC.out.report
        .map { report -> [ 'multiqc_done', report ] }

    //
    // MODULE: rf-count — per-base RT-stop or mutation counts from deduplicated BAM
    //
    ch_rfcount_with_fasta = ch_markdup_bam_bai
        .map { meta, bam, bai ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, bam, bai] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { _genome_build, bam_bai_tuple, fasta_tuple ->
            [ 'multiqc_done', [ bam_bai_tuple, fasta_tuple ] ]
        }
        .join(ch_multiqc_done)
        .map { _gate, rfcount_input, _multiqc_report ->
            rfcount_input
        }

    RNAFRAMEWORK_RFCOUNT (
        ch_rfcount_with_fasta.map { bam_bai_tuple, _fasta_tuple -> bam_bai_tuple },
        ch_rfcount_with_fasta.map { _bam_bai_tuple, fasta_tuple -> fasta_tuple }
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
            def group = "${meta.cell_line}__${meta.replicate}".toString()
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
    RNAFRAMEWORK_RFFOLD (
        RNAFRAMEWORK_RFNORM.out.xml
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
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    mapped_bam     = ch_dedup_bam                // channel: [ val(meta), path(bam) ]
    normalized_xml = RNAFRAMEWORK_RFNORM.out.xml // channel: [ val(meta), path(xml) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
