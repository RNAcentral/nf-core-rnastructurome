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
include { SAMTOOLS_FAIDX        } from '../modules/local/samtools/faidx/main'
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
    // MODULE: Run FastQC on raw reads (pre-trim)
    //
    FASTQC_PRE (
        ch_pretrim_fastqc_input
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    // Branch by probing principle so RT-stop and MaP can use different default cutadapt args
    def principle_branches = ch_samplesheet_for_branching.branch { meta, reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    def ch_reads_for_umi = principle_branches.rtstop.mix(principle_branches.map)
    def ch_reads_without_umi = ch_reads_for_umi.filter { meta, reads ->
        !((meta.umi_pattern ?: '').toString().trim())
    }

    UMITOOLS_EXTRACT (
        ch_reads_for_umi
    )

    def ch_reads_after_umi = ch_reads_without_umi.mix(UMITOOLS_EXTRACT.out.reads)
    def reads_for_cutadapt = ch_reads_after_umi.branch { meta, reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }
    def ch_rtstop_reads_for_cutadapt = reads_for_cutadapt.rtstop
    def ch_map_reads_for_cutadapt = reads_for_cutadapt.map
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.log.collect { umi_log -> umi_log[1] })

    //
    // MODULE: Run cutadapt
    //
    CUTADAPT_RTSTOP (
        ch_rtstop_reads_for_cutadapt
    )

    CUTADAPT_MAP (
        ch_map_reads_for_cutadapt
    )

    ch_trimmed_reads = CUTADAPT_RTSTOP.out.reads.mix(CUTADAPT_MAP.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_RTSTOP.out.log.collect { cutadapt_log -> cutadapt_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT_MAP.out.log.collect { cutadapt_log -> cutadapt_log[1] })

    //
    // MODULE: Run FastQC on trimmed reads (post-trim)
    //
    FASTQC_POST (
        ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    ch_rtstop_genome_build_fasta = principle_branches.rtstop
        .map { meta, reads ->
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
        .map { meta, reads ->
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
    ch_reference_fasta_keyed = ch_all_genome_build_fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }

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
    // MODULE: Index reference FASTA for coordinate-aware samtools sort
    //
    SAMTOOLS_FAIDX (
        ch_all_genome_build_fasta
    )

    ch_reference_fasta_fai_keyed = SAMTOOLS_FAIDX.out.fasta_fai
        .map { meta, fasta, fai -> [ meta.id.toString(), [meta, fasta, fai] ] }

    ch_sorted_inputs = ch_mapped_bam
        .map { meta, bam ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}' while preparing samtools sort input.")
            }
            [ genome_build, [meta, bam] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { genome_build, bam_tuple, fasta_tuple ->
            [ bam_tuple, fasta_tuple ]
        }

    SAMTOOLS_SORT (
        ch_sorted_inputs.map { bam_tuple, fasta_tuple -> bam_tuple },
        ch_sorted_inputs.map { bam_tuple, fasta_tuple -> fasta_tuple },
        false
    )

    SAMTOOLS_INDEX_SORT (
        SAMTOOLS_SORT.out.bam
    )

    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_SORT.out.bai.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
        .map { sample_id, bam_tuple, bai_tuple ->
            [ bam_tuple[0], bam_tuple[1], bai_tuple[1] ]
        }

    def dedup_branches = ch_sorted_bam_bai.branch { meta, bam, bai ->
        umi:     (meta.umi_pattern ?: '').toString().trim()
        non_umi: !((meta.umi_pattern ?: '').toString().trim())
    }

    UMITOOLS_DEDUP (
        dedup_branches.umi,
        false
    )

    ch_markdup_inputs = dedup_branches.non_umi
        .map { meta, bam, bai ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            if (!genome_build) {
                error("Missing genome_build for sample '${meta.id}' while preparing samtools markdup input.")
            }
            [ genome_build, [meta, bam] ]
        }
        .join(ch_reference_fasta_fai_keyed)
        .map { genome_build, bam_tuple, fasta_fai_tuple ->
            [ bam_tuple, fasta_fai_tuple ]
        }

    SAMTOOLS_MARKDUP (
        ch_markdup_inputs.map { bam_tuple, fasta_fai_tuple -> bam_tuple },
        ch_markdup_inputs.map { bam_tuple, fasta_fai_tuple -> fasta_fai_tuple }
    )

    ch_dedup_bam = UMITOOLS_DEDUP.out.bam.mix(SAMTOOLS_MARKDUP.out.bam)
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_DEDUP.out.log.collect { dedup_log -> dedup_log[1] })

    SAMTOOLS_INDEX_FINAL (
        ch_dedup_bam
    )

    ch_markdup_bam_bai = ch_dedup_bam
        .map { meta, bam -> [ meta.id.toString(), [meta, bam] ] }
        .join(SAMTOOLS_INDEX_FINAL.out.bai.map { meta, bai -> [ meta.id.toString(), [meta, bai] ] })
        .map { sample_id, bam_tuple, bai_tuple ->
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
        .map { genome_build, bam_bai_tuple, fasta_tuple ->
            [ bam_bai_tuple, fasta_tuple ]
        }

    SAMTOOLS_STATS (
        ch_stats_inputs.map { bam_bai_tuple, fasta_tuple -> bam_bai_tuple },
        ch_stats_inputs.map { bam_bai_tuple, fasta_tuple -> fasta_tuple }
    )

    SAMTOOLS_FLAGSTAT (
        ch_markdup_bam_bai
    )

    SAMTOOLS_IDXSTATS (
        ch_markdup_bam_bai
    )

    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { stats_file -> stats_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect { flagstat_file -> flagstat_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect { idxstats_file -> idxstats_file[1] })

    //
    // MODULE: rf-count — per-base RT-stop or mutation counts from deduplicated BAM
    //
    ch_rfcount_with_fasta = ch_markdup_bam_bai
        .map { meta, bam, bai ->
            def genome_build = (meta.genome_build ?: params.genome_build ?: params.genome)?.toString()
            [ genome_build, [meta, bam, bai] ]
        }
        .join(ch_reference_fasta_keyed)
        .map { genome_build, bam_bai_tuple, fasta_tuple ->
            [ bam_bai_tuple, fasta_tuple ]
        }

    RNAFRAMEWORK_RFCOUNT (
        ch_rfcount_with_fasta.map { bam_bai_tuple, fasta_tuple -> bam_bai_tuple },
        ch_rfcount_with_fasta.map { bam_bai_tuple, fasta_tuple -> fasta_tuple }
    )

    //
    // MODULE: rf-norm — normalise RC files to per-base reactivities (XML)
    //
    // Group RC files by 'group' key (defaults to sample_id) then join by condition.
    // TODO (rnacentral-probing-metadata-main): merge_metadata.py must populate
    // 'condition' (treated/untreated/denatured) and 'group' columns in the
    // generated samplesheet CSVs so that rf-norm can correctly pair samples.
    ch_rc_by_group = RNAFRAMEWORK_RFCOUNT.out.rc
        .map { meta, rc ->
            def group     = (meta.group ?: meta.sample_id ?: meta.id).toString()
            def condition = (meta.condition ?: 'treated').toLowerCase()
            [ group, condition, meta, rc ]
        }

    ch_treated   = ch_rc_by_group
        .filter  { group, condition, meta, rc -> condition == 'treated' }
        .map     { group, condition, meta, rc -> [ group, rc ] }
        .groupTuple()

    ch_untreated = ch_rc_by_group
        .filter  { group, condition, meta, rc -> condition == 'untreated' }
        .map     { group, condition, meta, rc -> [ group, rc ] }

    ch_denatured = ch_rc_by_group
        .filter  { group, condition, meta, rc -> condition == 'denatured' }
        .map     { group, condition, meta, rc -> [ group, rc ] }

    // Borrow principle/genome_build from the first sample in each group for group-level meta
    ch_group_meta = ch_rc_by_group
        .map     { group, condition, meta, rc -> [ group, meta ] }
        .groupTuple()
        .map     { group, metas -> [ group, metas[0] ] }

    ch_norm_input = ch_treated
        .join(ch_group_meta)
        .join(ch_untreated, remainder: true)
        .join(ch_denatured, remainder: true)
        .map { group, treated_rcs, base_meta, untreated_rc, denatured_rc ->
            def gmeta = base_meta + [ id: group ]
            [ gmeta, treated_rcs, untreated_rc ?: [], denatured_rc ?: [] ]
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
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    mapped_bam     = ch_dedup_bam                // channel: [ val(meta), path(bam) ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
