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
include { ENSEMBL_GTF          } from '../modules/local/ensembl/gtf/main'
include { NCBI_FASTA           } from '../modules/local/ncbi/fasta/main'
include { NCBI_GTF             } from '../modules/local/ncbi/gtf/main'
include { FASTA_SORT as FASTA_SORT_LOCAL    } from '../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_ENSEMBL } from '../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_NCBI    } from '../modules/local/fasta/sort/main'
include { RNAFRAMEWORK_DOTPLOT2BP } from '../modules/local/dotplot2bp/main'
include { MERGE_BP               } from '../modules/local/merge_bp/main'
include { RNAFRAMEWORK_RFWIGGLE  } from '../modules/local/rnaframework/wiggle/main'
include { MERGE_WIG              } from '../modules/local/merge_wig/main'
include { AVERAGE_WIG            } from '../modules/local/average_wig/main'
include { MERGE_SHANNON_WIG      } from '../modules/local/merge_shannon_wig/main'
include { RNAFRAMEWORK_TORDAT    } from '../modules/local/tordat/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_REACTIVITY } from '../modules/nf-core/ucsc/wigtobigwig/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_SHANNON    } from '../modules/nf-core/ucsc/wigtobigwig/main'
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
    ch_samplesheet         // channel: samplesheet read in from --input
    pipeline_config_input  // map: pipeline configuration captured at the entry workflow
    main:

    def pipeline_config = defaultPipelineConfig() + (pipeline_config_input ?: [:])
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

    def ch_pretrim_reads_split = CAT_FASTQ.out.reads.multiMap { meta, reads ->
        fastqc: [ meta, reads ]
        branching: [ meta, reads ]
    }
    def ch_pretrim_fastqc_input = ch_pretrim_reads_split.fastqc.map { meta, reads -> [ meta, reads ] }
    def ch_samplesheet_for_branching = ch_pretrim_reads_split.branching.map { meta, reads -> [ meta, reads ] }

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
    def ch_reads_with_umi = ch_reads_for_umi.filter { meta, _reads ->
        (meta.umi_pattern ?: '').toString().trim()
    }
    def ch_reads_without_umi = ch_reads_for_umi.filter { meta, _reads ->
        !((meta.umi_pattern ?: '').toString().trim())
    }

    //
    // MODULE: umi_tools extract — extract UMIs from reads
    //
    UMITOOLS_EXTRACT (
        ch_reads_with_umi
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

    def ch_rtstop_trimmed_split = CUTADAPT_RTSTOP.out.reads.multiMap { meta, reads ->
        align: [ meta, reads ]
        fastqc: [ meta, reads ]
    }
    def ch_map_trimmed_split = CUTADAPT_MAP.out.reads.multiMap { meta, reads ->
        align: [ meta, reads ]
        fastqc: [ meta, reads ]
    }
    def ch_rtstop_trimmed_for_align = ch_rtstop_trimmed_split.align.map { meta, reads -> [ meta, reads ] }
    def ch_rtstop_trimmed_for_fastqc = ch_rtstop_trimmed_split.fastqc.map { meta, reads -> [ meta, reads ] }
    def ch_map_trimmed_for_align = ch_map_trimmed_split.align.map { meta, reads -> [ meta, reads ] }
    def ch_map_trimmed_for_fastqc = ch_map_trimmed_split.fastqc.map { meta, reads -> [ meta, reads ] }
    ch_trimmed_reads = ch_rtstop_trimmed_for_fastqc.mix(ch_map_trimmed_for_fastqc)
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
        ch_trimmed_reads.map { meta, reads -> [ meta + [id: "${meta.id}_trimmed"], reads ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect { fastqc_zip -> fastqc_zip[1] })

    ch_reference_requests = ch_samplesheet_for_branching
        .map { meta, _reads ->
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def original_organism = meta.organism?.toString() ?: reference_key
            if (pipeline_config.fasta) {
                return [ reference_key, "path::${pipeline_config.fasta.toString()}", original_organism ]
            }

            def genome_entry = pipeline_config.genomes?.containsKey(reference_key) ? pipeline_config.genomes[reference_key] : null
            def transcript_fasta = genome_entry?.transcript_fasta ?: genome_entry?.transcriptome ?: genome_entry?.cdna
            if (transcript_fasta) {
                return [ reference_key, "path::${transcript_fasta.toString()}", original_organism ]
            }

            def explicit_ensembl = genome_entry?.ensembl_species ?: pipeline_config.ensembl_species_map?.get(reference_key)
            if (explicit_ensembl) {
                return [ reference_key, "ensembl::${explicit_ensembl.toLowerCase()}", original_organism ]
            }

            // Pre-configured NCBI accessions (e.g. from viral_genomes.config) skip the
            // Ensembl step entirely — faster and avoids spurious 404 log lines.
            def ncbi_accessions = genome_entry?.ncbi_accessions ?: pipeline_config.ncbi_accessions_map?.get(reference_key)
            if (ncbi_accessions) {
                def acc_str = (ncbi_accessions instanceof List) ? ncbi_accessions.join(',') : ncbi_accessions.toString()
                return [ reference_key, "ncbi::${acc_str}", original_organism ]
            }

            // All remaining organisms go to Ensembl first.  If Ensembl returns HTTP 404 the
            // module emits a not_found signal which is automatically routed to NCBI_FASTA.
            if (!reference_key) {
                error("No organism specified for sample '${meta.id}'. Provide --fasta, --organism, or set params.genomes.")
            }
            return [ reference_key, "ensembl::${reference_key}", original_organism ]
        }
        .groupTuple()
        .map { reference_key, resolutions, organisms ->
            def unique_resolutions = resolutions.unique()
            if (unique_resolutions.size() != 1) {
                error("Multiple transcript reference resolutions were detected for reference '${reference_key}': ${unique_resolutions.join(', ')}")
            }
            [ reference_key, unique_resolutions[0], organisms[0] ]
        }

    ch_reference_local = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('path::') }
        .map { reference_key, resolution, _original_organism ->
            def fasta_path = resolution - 'path::'
            [ [ id: reference_key, organism: reference_key ], file(fasta_path, checkIfExists: true) ]
        }

    ch_reference_ensembl_input = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ensembl::') }
        .map { reference_key, resolution, original_organism ->
            def ensembl_species = resolution - 'ensembl::'
            [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species,
                original_organism: original_organism ], ensembl_species ]
        }

    ENSEMBL_TRANSCRIPTOME (
        ch_reference_ensembl_input,
        [
            ensembl_release : pipeline_config.ensembl_release,
            ensembl_base_url: pipeline_config.ensembl_base_url
        ],
        file("${projectDir}/bin/ensembl_transcriptome.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(ENSEMBL_TRANSCRIPTOME.out.versions)

    // Organisms not found on Ensembl FTP are routed to NCBI_FASTA for automatic accession
    // search.  meta.original_organism carries the raw samplesheet organism string used as
    // the esearch query when no accessions are pre-configured.
    ch_reference_ncbi_from_ensembl = ENSEMBL_TRANSCRIPTOME.out.not_found
        .map { meta, _not_found_file -> [ meta, "" ] }

    ch_reference_ncbi_explicit = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ncbi::') }
        .map { reference_key, resolution, original_organism ->
            def accessions = resolution - 'ncbi::'
            [ [ id: reference_key, organism: reference_key, original_organism: original_organism ], accessions ]
        }

    ch_reference_ncbi_input = ch_reference_ncbi_explicit.mix(ch_reference_ncbi_from_ensembl)

    NCBI_FASTA (
        ch_reference_ncbi_input,
        file("${projectDir}/bin/ncbi_fasta.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(NCBI_FASTA.out.versions)

    NCBI_GTF (
        NCBI_FASTA.out.fasta,
        file("${projectDir}/bin/ncbi_gtf.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(NCBI_GTF.out.versions)

    ch_reference_gtf_requests = ch_samplesheet_for_branching
        .map { meta, _reads ->
            def reference_key    = resolveReferenceKey(meta, pipeline_config.organism)
            def original_organism = meta.organism?.toString() ?: reference_key

            if (pipeline_config.gtf) {
                return [ reference_key, "path::${pipeline_config.gtf.toString()}", original_organism ]
            }

            def genome_entry = pipeline_config.genomes?.containsKey(reference_key) ? pipeline_config.genomes[reference_key] : null
            def gtf_path = genome_entry?.gtf
            if (gtf_path) {
                return [ reference_key, "path::${gtf_path.toString()}", original_organism ]
            }

            def explicit_ensembl = genome_entry?.ensembl_species ?: pipeline_config.ensembl_species_map?.get(reference_key)
            if (explicit_ensembl) {
                return [ reference_key, "ensembl::${explicit_ensembl.toLowerCase()}", original_organism ]
            }

            // NCBI references (pre-configured or auto-search via Ensembl 404 fallback) do not
            // go through the Ensembl GTF route.  Their annotation is generated by NCBI_GTF,
            // which derives a synthetic transcript record per accession from NCBI_FASTA output.
            def ncbi_accessions = genome_entry?.ncbi_accessions ?: pipeline_config.ncbi_accessions_map?.get(reference_key)
            if (ncbi_accessions) {
                return [ reference_key, "none::", original_organism ]
            }

            return [ reference_key, "ensembl::${reference_key}", original_organism ]
        }
        .groupTuple()
        .map { reference_key, resolutions, organisms ->
            def unique_resolutions = resolutions.unique()
            if (unique_resolutions.size() != 1) {
                error("Multiple GTF reference resolutions were detected for reference '${reference_key}': ${unique_resolutions.join(', ')}")
            }
            [ reference_key, unique_resolutions[0], organisms[0] ]
        }

    // Local GTF (already in GTF format) — used directly
    ch_reference_gtf_local = ch_reference_gtf_requests
        .filter { _k, resolution, _o -> resolution.startsWith('path::') }
        .map { reference_key, resolution, _original_organism ->
            [ [ id: reference_key, organism: reference_key ], file(resolution - 'path::', checkIfExists: true) ]
        }

    ch_reference_gtf_ensembl_input = ch_reference_gtf_requests
        .filter { _k, resolution, _o -> resolution.startsWith('ensembl::') }
        .map { reference_key, resolution, original_organism ->
            def ensembl_species = resolution - 'ensembl::'
            [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species,
                original_organism: original_organism ], ensembl_species ]
        }

    ENSEMBL_GTF (
        ch_reference_gtf_ensembl_input,
        [
            ensembl_release : pipeline_config.ensembl_release,
            ensembl_base_url: pipeline_config.ensembl_base_url
        ],
        file("${projectDir}/bin/ensembl_gtf.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(ENSEMBL_GTF.out.versions)

    // NCBI references (both pre-configured and Ensembl-not-found): annotation is the synthetic
    // GTF from NCBI_GTF, which has already run above from NCBI_FASTA.out.fasta.
    // ENSEMBL_GTF.out.not_found is silently dropped — those organisms have a GTF via NCBI_GTF.
    ch_all_reference_gtf = ch_reference_gtf_local
        .mix(ENSEMBL_GTF.out.gtf)
        .mix(NCBI_GTF.out.gtf)

    FASTA_SORT_LOCAL (
        ch_reference_local
    )
    ch_versions = ch_versions.mix(FASTA_SORT_LOCAL.out.versions)

    FASTA_SORT_ENSEMBL (
        ENSEMBL_TRANSCRIPTOME.out.fasta
    )
    ch_versions = ch_versions.mix(FASTA_SORT_ENSEMBL.out.versions)

    FASTA_SORT_NCBI (
        NCBI_FASTA.out.fasta
    )
    ch_versions = ch_versions.mix(FASTA_SORT_NCBI.out.versions)

    ch_reference_fasta_keyed = FASTA_SORT_LOCAL.out.fasta
        .mix(FASTA_SORT_ENSEMBL.out.fasta)
        .mix(FASTA_SORT_NCBI.out.fasta)
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    ch_reference_gtf_keyed   = ch_all_reference_gtf.map { meta, gtf -> [ meta.id.toString(), [meta, gtf] ] }
    ch_reference_fasta_map   = ch_reference_fasta_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
    ch_reference_gtf_map     = ch_reference_gtf_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }

    ch_rtstop_reference_fasta = principle_branches.rtstop
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def ref_map = combined[2]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def ref_tuple = ref_map[reference_key]
            if (!ref_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' in RT-stop branch.")
            }
            ref_tuple
        }
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        .groupTuple()
        .map { _reference_key, entries -> entries[0] }

    ch_map_reference_fasta = principle_branches.map
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def ref_map = combined[2]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def ref_tuple = ref_map[reference_key]
            if (!ref_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' in MaP branch.")
            }
            ref_tuple
        }
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        .groupTuple()
        .map { _reference_key, entries -> entries[0] }

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
    ch_bowtie_index_map = ch_bowtie_index_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
    ch_bowtie2_index_map = ch_bowtie2_index_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }

    ch_rtstop_align_inputs = ch_rtstop_trimmed_for_align
        .combine(ch_bowtie_index_map)
        .map { combined ->
            def meta = combined[0]
            def reads = combined[1]
            def index_map = combined[2]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def index_tuple = index_map[reference_key]
            if (!index_tuple) {
                error("No Bowtie index resolved for reference '${reference_key}'.")
            }
            [ [meta, reads], index_tuple ]
        }

    ch_map_align_inputs = ch_map_trimmed_for_align
        .combine(ch_bowtie2_index_map)
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def reads = combined[1]
            def index_map = combined[2]
            def ref_map = combined[3]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def index_tuple = index_map[reference_key]
            def fasta_tuple = ref_map[reference_key]
            if (!index_tuple) {
                error("No Bowtie2 index resolved for reference '${reference_key}'.")
            }
            if (!fasta_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' in MaP alignment.")
            }
            [ [meta, reads], index_tuple, fasta_tuple ]
        }

    //
    // MODULE: bowtie align — align RT-stop reads with Bowtie v1
    //
    def ch_rtstop_align_split = ch_rtstop_align_inputs.multiMap { entry ->
        reads: entry[0]
        index: entry[1]
    }
    def ch_rtstop_align_reads = ch_rtstop_align_split.reads
    def ch_rtstop_align_index = ch_rtstop_align_split.index
    BOWTIE_ALIGN (
        ch_rtstop_align_reads,
        ch_rtstop_align_index,
        false
    )

    //
    // MODULE: bowtie2 align — align MaP reads with Bowtie2
    //
    def ch_map_align_split = ch_map_align_inputs.multiMap { entry ->
        reads: entry[0]
        index: entry[1]
        fasta: entry[2]
    }
    def ch_map_align_reads = ch_map_align_split.reads
    def ch_map_align_index = ch_map_align_split.index
    def ch_map_align_fasta = ch_map_align_split.fasta
    BOWTIE2_ALIGN (
        ch_map_align_reads,
        ch_map_align_index,
        ch_map_align_fasta,
        false,
        false
    )

    ch_mapped_bam = BOWTIE_ALIGN.out.bam.mix(BOWTIE2_ALIGN.out.bam)
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect { bowtie_log -> bowtie_log[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { bowtie2_log -> bowtie2_log[1] })

    //
    // MODULE: samtools sort — coordinate-sort mapped BAMs
    // FASTA/FAI not needed for BAM output (only required for CRAM); pass empty.
    //
    SAMTOOLS_SORT (
        ch_mapped_bam,
        channel.value([ [], [], [] ]),
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
    // MODULE: samtools faidx — index reference FASTA for markdup
    //
    SAMTOOLS_FAIDX (
        FASTA_SORT_LOCAL.out.fasta.mix(FASTA_SORT_ENSEMBL.out.fasta).map { meta, fasta -> [meta, fasta, []] },
        false
    )

    //
    // MODULE: samtools markdup — deduplicate non-UMI BAMs
    // FASTA/FAI not needed for BAM output; pass empty.
    // Drop the bai from [meta, bam, bai] — markdup only takes [meta, bam].
    //
    SAMTOOLS_MARKDUP (
        dedup_branches.non_umi.map { meta, bam, _bai -> [ meta, bam ] },
        channel.value([ [], [], [] ])
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
    // MODULE: rf-count — per-base RT-stop or mutation counts from deduplicated BAM
    //
    ch_rfcount_with_fasta = ch_markdup_bam_bai
        .combine(ch_reference_fasta_map)
        .map { combined ->
            def meta = combined[0]
            def bam = combined[1]
            def bai = combined[2]
            def ref_map = combined[3]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def fasta_tuple = ref_map[reference_key]
            if (!fasta_tuple) {
                error("No transcript FASTA resolved for reference '${reference_key}' for rf-count.")
            }
            [ [meta, bam, bai], fasta_tuple ]
        }

    def ch_rfcount_split = ch_rfcount_with_fasta.multiMap { entry ->
        bam: entry[0]
        fasta: entry[1]
    }
    def ch_rfcount_bam_input = ch_rfcount_split.bam
    def ch_rfcount_fasta_input = ch_rfcount_split.fasta
    RNAFRAMEWORK_RFCOUNT (
        ch_rfcount_bam_input,
        ch_rfcount_fasta_input
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

    RNAFRAMEWORK_RFFOLD (
        ch_fold_input
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

    //
    // MODULE: merge_bp — merge per-transcript .bp files into a single file per fold group for genome browser visualisation
    //
    MERGE_BP (
        RNAFRAMEWORK_DOTPLOT2BP.out.bp
    )

    //
    // MODULE: rf-wiggle — convert rf-norm XML reactivities to WIG + chrom.sizes
    //
    RNAFRAMEWORK_RFWIGGLE (
        RNAFRAMEWORK_RFNORM.out.xml
    )

    def ch_merge_wiggle_input = RNAFRAMEWORK_RFWIGGLE.out.wig
        .map { meta, wig -> [ meta.id.toString(), [meta, wig] ] }
        .join(RNAFRAMEWORK_RFNORM.out.xml.map { meta, xml -> [ meta.id.toString(), [meta, xml] ] })
        .map { _sample_id, wig_tuple, xml_tuple ->
            [ wig_tuple[0], wig_tuple[1], xml_tuple[1] ]
        }

    MERGE_WIG (
        ch_merge_wiggle_input
    )

    //
    // Group per-replicate merged WIGs by cell_line.
    // Single replicate: bypass AVERAGE_WIG, keep meta.id = "HEK293T_1" → HEK293T_1_reactivity.bw
    // Multiple replicates: run AVERAGE_WIG, set meta.id = cell_line   → HEK293T_reactivity.bw
    //
    def ch_reactivity_grouped = MERGE_WIG.out.merged_wig
        .map { meta, wig -> [ meta.id.toString(), meta, wig ] }
        .join(MERGE_WIG.out.chrom_sizes.map { meta, sizes -> [ meta.id.toString(), sizes ] })
        .map { _id, meta, wig, sizes -> [ meta.cell_line.toString(), meta, wig, sizes ] }
        .groupTuple(by: 0)
        .map { cell_line, metas, wigs, sizes_list ->
            def base_meta = wigs.size() == 1
                ? metas[0]
                : metas[0] + [ id: cell_line ]
            [ base_meta, wigs.flatten(), sizes_list[0] ]
        }

    def ch_reactivity_branches = ch_reactivity_grouped.branch { _meta, wigs, _sizes ->
        single: wigs.size() == 1
        multi:  true
    }

    def ch_single_reactivity = ch_reactivity_branches.single.multiMap { meta, wigs, sizes ->
        wig:   [ meta, wigs[0] ]
        sizes: sizes
    }

    AVERAGE_WIG (
        ch_reactivity_branches.multi
    )

    //
    // MODULE: wigToBigWig — convert merged/averaged WIG to BigWig for IGV
    //
    UCSC_WIGTOBIGWIG_REACTIVITY (
        ch_single_reactivity.wig.mix(AVERAGE_WIG.out.merged_wig),
        ch_single_reactivity.sizes.mix(AVERAGE_WIG.out.chrom_sizes.map { _meta, sizes -> sizes })
    )

    //
    // MODULE: merge_shannon_wig + wigToBigWig — merge per-transcript Shannon entropy WIG files
    // and convert to BigWig for genome browser visualisation
    //
    if (params.rffold_shannon_entropy) {
        def ch_shannon_merge_input = RNAFRAMEWORK_RFFOLD.out.shannon_wig
            .map { meta, wigs -> [ meta.id.toString(), meta, wigs ] }
            .join(ch_fold_input.map { meta, xmls -> [ meta.id.toString(), xmls ] })
            .map { _id, meta, wigs, xmls -> [ meta, wigs, xmls ] }

        MERGE_SHANNON_WIG (
            ch_shannon_merge_input
        )

        UCSC_WIGTOBIGWIG_SHANNON (
            MERGE_SHANNON_WIG.out.merged_wig,
            MERGE_SHANNON_WIG.out.chrom_sizes.map { _meta, sizes -> sizes }
        )
    }

    //
    // MODULE: tordat — compile rf-norm XML + rf-fold .db structures into RDAT format
    //
    def ch_rdat_input = ch_fold_input
        .map { meta, xml -> [ meta.id.toString(), meta, xml ] }
        .combine(
            RNAFRAMEWORK_RFFOLD.out.structures
                .map { meta, fold_dir -> [ meta.id.toString(), fold_dir ] },
            by: 0
        )
        .map { _key, fold_meta, xml, fold_dir -> [ fold_meta, xml, fold_dir ] }
        .map { fold_meta, xml, fold_dir ->
            def reference_key = resolveReferenceKey(fold_meta, pipeline_config.organism)
            def fasta_name = pipeline_config.fasta
                ? file(pipeline_config.fasta.toString()).name
                : "${reference_key}.transcripts.fa.gz"
            [ fold_meta + [ fasta_name: fasta_name ], xml, fold_dir ]
        }

    RNAFRAMEWORK_TORDAT (
        ch_rdat_input,
        file("${projectDir}/bin/rnaframework_to_rdat.py", checkIfExists: true)
    )

    // Add RNAframework outputs to MultiQC input collection.
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFCOUNT.out.rc.collect { rc_file -> rc_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFCOUNT.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.xml.collect { xml_file -> xml_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFNORM.out.plots.collect { plot_file -> plot_file[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RNAFRAMEWORK_RFFOLD.out.structures.collect { fold_dir -> fold_dir[1] })

    // RF-norm summary table: one row per normalisation group (cell_line + replicate).
    def ch_rfnorm_stats_mqc = RNAFRAMEWORK_RFNORM.out.log
        .map { meta, log -> [ meta.id.toString(), parseRfnormLog(log) ] }
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
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFCOUNT.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFFOLD.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_DOTPLOT2BP.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BP.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFWIGGLE.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_WIG.out.versions.first())
    ch_versions = ch_versions.mix(AVERAGE_WIG.out.versions.first())
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
    multiqc_report   = MULTIQC.out.report.toList()        // channel: /path/to/multiqc_report.html
    mapped_bam       = ch_dedup_bam                       // channel: [ val(meta), path(bam) ]
    normalized_xml   = RNAFRAMEWORK_RFNORM.out.xml        // channel: [ val(meta), path(xml) ]
    fold_structures  = RNAFRAMEWORK_RFFOLD.out.structures // channel: [ val(meta), path(dir) ]
    fold_bp          = RNAFRAMEWORK_DOTPLOT2BP.out.bp     // channel: [ val(meta), path(bp) ]
    merged_bp        = MERGE_BP.out.bp                   // channel: [ val(meta), path(*_merged.bp) ]
    versions         = ch_versions                        // channel: [ path(versions.yml) ]

}

def normaliseEnsemblSpecies(value) {
    value
        ?.toString()
        ?.trim()
        ?.toLowerCase()
        ?.replaceAll(/[^a-z0-9_]+/, '_')
        ?.replaceAll(/^_+|_+$/, '')
}

def defaultPipelineConfig() {
    [
        organism                          : null,
        fasta                             : null,
        gtf                               : null,
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
        rfnorm_nan                        : 10,
        rnaframework_r_path               : '/usr/bin/R'
    ]
}

def resolveReferenceKey(meta, fallbackOrganism) {
    def sampleId = meta?.id ?: 'unknown'
    def rawReference = (meta?.organism ?: fallbackOrganism)?.toString()?.trim()
    if (!rawReference) {
        error("Missing organism for sample '${sampleId}'. Set organism in the samplesheet or provide --organism.")
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
        'ncbi_accessions_map',
        'genomes',
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

        def filteredSection = sectionParams.findAll { key, _value ->
            if (hiddenKeys.contains(key)) {
                return false
            }
            true
        }

        [(sectionName): filteredSection]
    }.findAll { _sectionName, sectionParams ->
        !(sectionParams instanceof Map) || !sectionParams.isEmpty()
    }
}

def addModuleOptionsSummary(summaryParams, pipeline_config) {
    def sampleMetadata = parseInputSamplesheetMetadata(pipeline_config.input)
    def moduleOptions = buildModuleOptionsSummary(pipeline_config, sampleMetadata)
    if (moduleOptions.isEmpty()) {
        return summaryParams
    }
    summaryParams + ['Module options': moduleOptions]
}

def parseInputSamplesheetMetadata(inputPath) {
    if (!inputPath) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def inputFile = file(inputPath.toString())
    if (!inputFile.exists()) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def lines = inputFile.readLines().findAll { line -> line?.trim() }
    if (lines.size() < 2) {
        return [principles: [], conditions: [], methods: [], adapter_5p: [], adapter_3p: []]
    }

    def header = lines[0].split(',', -1)*.trim()
    def principleIdx = header.indexOf('principle')
    def conditionIdx = header.indexOf('condition')
    def methodIdx = header.indexOf('method')
    def adapter5pIdx = header.indexOf('adapter_5p')
    def adapter3pIdx = header.indexOf('adapter_3p')

    def principles = []
    def conditions = []
    def methods = []
    def adapter5p = []
    def adapter3p = []

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
        if (adapter5pIdx >= 0 && adapter5pIdx < fields.size()) {
            def value = fields[adapter5pIdx]?.trim()
            if (value) adapter5p << value
        }
        if (adapter3pIdx >= 0 && adapter3pIdx < fields.size()) {
            def value = fields[adapter3pIdx]?.trim()
            if (value) adapter3p << value
        }
    }

    [
        principles: principles.unique(),
        conditions: conditions.collect { condition -> condition.toLowerCase() }.unique(),
        methods   : methods.unique(),
        adapter_5p: adapter5p.unique(),
        adapter_3p: adapter3p.unique()
    ]
}

def buildModuleOptionsSummary(pipeline_config, sampleMetadata) {
    def moduleOptions = [:]
    def principles = (sampleMetadata.principles ?: []).collect { principle -> principle.toLowerCase() }
    def adapter5p = (sampleMetadata.adapter_5p ?: []).findAll { adapter -> adapter?.trim() }
    def adapter3p = (sampleMetadata.adapter_3p ?: []).findAll { adapter -> adapter?.trim() }

    if (!adapter5p.isEmpty()) {
        moduleOptions['cutadapt_adapter_5p'] = adapter5p.join(', ')
    }
    if (!adapter3p.isEmpty()) {
        moduleOptions['cutadapt_adapter_3p'] = adapter3p.join(', ')
    }

    if (!principles || principles.contains('rt-stop')) {
        moduleOptions['bowtie_rtstop_aligner'] = 'bowtie'
        moduleOptions['bowtie_rtstop_args'] = renderBowtie1Args(pipeline_config)
    }
    if (principles.contains('map')) {
        moduleOptions['bowtie_map_aligner'] = 'bowtie2'
        moduleOptions['bowtie_map_args'] = renderBowtie2Args(pipeline_config)
    }

    def rfnormSummary = renderRfNormSummary(pipeline_config, sampleMetadata)
    moduleOptions.putAll(rfnormSummary)

    moduleOptions.findAll { _k, v -> v != null && v.toString().trim() }
}

def renderBowtie1Args(pipeline_config) {
    def manualOnly = pipeline_config.bowtie_manual_only as Boolean ?: false
    def manualParams = (pipeline_config.bowtie_mapping_params ?: '').toString().trim()
    if (manualOnly) {
        return manualParams ?: 'none'
    }
    def args = []
    if (pipeline_config.bowtie_all as Boolean) {
        args << '-a'
    } else if (pipeline_config.bowtie_k != null) {
        args << "-k ${pipeline_config.bowtie_k as Integer}"
    }
    if (pipeline_config.bowtie_norc as Boolean) {
        args << '--norc'
    }
    if ((pipeline_config.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${pipeline_config.bowtie_trim5 as Integer}"
    }
    if ((pipeline_config.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${pipeline_config.bowtie_trim3 as Integer}"
    }
    def seedlen = pipeline_config.bowtie_seedlen != null ? pipeline_config.bowtie_seedlen as Integer : 28
    args << "-l ${seedlen}"
    if (pipeline_config.bowtie_v != null) {
        args << "-v ${pipeline_config.bowtie_v as Integer}"
    } else {
        args << "-n ${pipeline_config.bowtie_n as Integer}"
    }
    if (!(pipeline_config.bowtie_all as Boolean) && pipeline_config.bowtie_k == null && pipeline_config.bowtie_max != null) {
        args << "-m ${pipeline_config.bowtie_max as Integer}"
    }
    args << "--chunkmbs ${pipeline_config.bowtie_chunkmbs as Integer}"
    if (manualParams) {
        args << manualParams
    }
    args.join(' ').trim()
}

def renderBowtie2Args(pipeline_config) {
    def manualOnly = pipeline_config.bowtie_manual_only as Boolean ?: false
    def manualParams = (pipeline_config.bowtie_mapping_params ?: '').toString().trim()
    if (manualOnly) {
        return manualParams ?: 'none'
    }
    def args = []
    if (pipeline_config.bowtie_all as Boolean) {
        args << '-a'
    } else if (pipeline_config.bowtie_k != null) {
        args << "-k ${pipeline_config.bowtie_k as Integer}"
    }
    if (pipeline_config.bowtie_norc as Boolean) {
        args << '--norc'
    }
    if ((pipeline_config.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${pipeline_config.bowtie_trim5 as Integer}"
    }
    if ((pipeline_config.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${pipeline_config.bowtie_trim3 as Integer}"
    }
    def seedlen = pipeline_config.bowtie_seedlen != null ? pipeline_config.bowtie_seedlen as Integer : 22
    args << "-L ${seedlen}"
    args << "-N ${pipeline_config.bowtie2_N as Integer}"
    args << "-D ${pipeline_config.bowtie2_D as Integer}"
    args << "-R ${pipeline_config.bowtie2_R as Integer}"
    args << "--mp ${pipeline_config.bowtie2_mp}"
    args << "--dpad ${pipeline_config.bowtie2_dpad as Integer}"
    args << "--rdg ${pipeline_config.bowtie2_rdg}"
    args << "--rfg ${pipeline_config.bowtie2_rfg}"
    if (pipeline_config.bowtie2_softclip as Boolean) {
        args << '--local'
        args << "--ma ${pipeline_config.bowtie2_ma as Integer}"
    }
    if (pipeline_config.bowtie2_dovetail as Boolean) {
        args << '--dovetail'
    }
    if (manualParams) {
        args << manualParams
    }
    args.join(' ').trim()
}

def renderRfNormSummary(pipeline_config, sampleMetadata) {
    def principles = (sampleMetadata.principles ?: []).collect { principle -> principle.toLowerCase() }.unique()
    def conditions = (sampleMetadata.conditions ?: []).collect { condition -> condition.toLowerCase() }.unique()
    if (principles.size() != 1) {
        return [
            rfnorm_mode: 'dynamic (mixed principles across samples)'
        ]
    }

    def principle = principles[0]
    def hasUntreated = conditions.contains('untreated')
    def hasDenatured = conditions.contains('denatured')
    def scoringMethod = principle == 'map' ? (hasUntreated ? 3 : 4) : (hasUntreated ? 1 : 2)
    def normMethod = resolveRfNormNormMethod(pipeline_config, scoringMethod)
    def isDmsOnly = ((sampleMetadata.methods ?: []).collect { method -> method.toLowerCase() }.unique()) == ['dms']
    def isDmsBroad = isDmsOnly && sampleMetadata.pH != null && (sampleMetadata.pH as Double) >= 8.0
    def reactiveBases = pipeline_config.rfnorm_reactive_bases ?: (isDmsOnly ? (isDmsBroad ? 'ACGU' : 'AC') : null)
    def dynamicWindow = pipeline_config.rfnorm_dynamic_window != null ? (pipeline_config.rfnorm_dynamic_window as Integer) : (isDmsOnly && !isDmsBroad ? 50 : null)

    def args = [
        "-sm ${scoringMethod}",
        "-nm ${normMethod}"
    ]
    if (pipeline_config.rfnorm_remap_reactivities as Boolean) args << '--remap-reactivities'
    if (reactiveBases) args << "--reactive-bases ${reactiveBases}"
    if (pipeline_config.rfnorm_norm_window != null) args << "--norm-window ${pipeline_config.rfnorm_norm_window as Integer}"
    if (pipeline_config.rfnorm_window_offset != null) args << "--window-offset ${pipeline_config.rfnorm_window_offset as Integer}"
    if (dynamicWindow != null) args << "--dynamic-window ${dynamicWindow}"
    if (pipeline_config.rfnorm_norm_independent as Boolean) args << '--norm-independent'
    if (pipeline_config.rfnorm_norm_factor) args << "--norm-factor ${pipeline_config.rfnorm_norm_factor}"
    if (pipeline_config.rfnorm_raw as Boolean) args << '--raw'
    if (pipeline_config.rfnorm_pseudocount != null) args << "--pseudocount ${pipeline_config.rfnorm_pseudocount}"
    if (pipeline_config.rfnorm_max_score != null) args << "--max-score ${pipeline_config.rfnorm_max_score}"
    if (pipeline_config.rfnorm_ignore_lower_than_untreated as Boolean) args << '--ignore-lower-than-untreated'
    if (pipeline_config.rfnorm_max_untreated_mut != null) args << "--max-untreated-mut ${pipeline_config.rfnorm_max_untreated_mut}"
    if (pipeline_config.rfnorm_max_mutation_rate != null) args << "--max-mutation-rate ${pipeline_config.rfnorm_max_mutation_rate}"
    def meanCoverage = pipeline_config.rfnorm_mean_coverage != null ? pipeline_config.rfnorm_mean_coverage as BigDecimal : 0
    if (meanCoverage > 0) args << "--mean-coverage ${pipeline_config.rfnorm_mean_coverage}"
    def medianCoverage = pipeline_config.rfnorm_median_coverage != null ? pipeline_config.rfnorm_median_coverage as BigDecimal : 0
    if (medianCoverage > 0) args << "--median-coverage ${pipeline_config.rfnorm_median_coverage}"
    def nanThreshold = pipeline_config.rfnorm_nan != null ? pipeline_config.rfnorm_nan as Integer : 10
    if (nanThreshold != 10) args << "--nan ${nanThreshold}"
    args << '--img'
    args << "-R ${pipeline_config.rnaframework_r_path}"

    [
        rfnorm_mode         : "${principle.toUpperCase()} ${hasUntreated ? 'with untreated' : 'treated-only'}${hasDenatured ? ' + denatured' : ''}",
        rfnorm_scoring      : "${rfNormScoringLabel(scoringMethod)} (sm=${scoringMethod})",
        rfnorm_normalisation: "${rfNormNormLabel(normMethod)} (nm=${normMethod})",
        rfnorm_args         : args.join(' ').trim()
    ]
}

def resolveRfNormNormMethod(pipeline_config, scoringMethod) {
    def defaultMethod = (scoringMethod as Integer) == 2 ? 2 : 3
    if (pipeline_config.rfnorm_norm_method == null) {
        return defaultMethod
    }

    def requestedMethod = pipeline_config.rfnorm_norm_method as Integer
    if (!(requestedMethod in [2, 3])) {
        error("Unsupported rf-norm normalization method '${pipeline_config.rfnorm_norm_method}'. Expected one of: 2, 3.")
    }

    requestedMethod
}

def rfNormScoringLabel(code) {
    def labels = [
        1: 'Ding',
        2: 'Rouskin',
        3: 'Siegfried',
        4: 'Zubradt'
    ]
    labels[code as Integer] ?: 'unknown'
}

def rfNormNormLabel(code) {
    def labels = [
        2: '90% Winsorizing',
        3: 'Box-plot normalisation'
    ]
    labels[code as Integer] ?: 'unknown'
}

def parseCutadaptCommandArg(logFile, optionName) {
    def commandLine = logFile.readLines().find { line -> line.startsWith('Command line parameters:') }
    if (!commandLine) {
        return 'none'
    }
    def pattern = java.util.regex.Pattern.compile("(?:^|\\s)${java.util.regex.Pattern.quote(optionName)}\\s+(\\S+)")
    def matcher = pattern.matcher(commandLine)
    matcher.find() ? matcher.group(1) : 'none'
}

def countProgressionMultiqc(rows) {
    def rowEntries
    if (rows instanceof Map) {
        rowEntries = rows.entrySet().collect { entry -> [entry.key, entry.value] }
    } else if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else if (rows instanceof List && rows.every { row -> row instanceof List && row.size() == 2 && row[1] instanceof Map }) {
        rowEntries = rows
    } else if (rows instanceof List && rows.size() % 2 == 0 && rows.collate(2).every { pair -> pair.size() == 2 && pair[1] instanceof Map }) {
        rowEntries = rows.collate(2)
    } else {
        error("Unexpected count progression row structure: ${rows?.getClass()?.name} -> ${rows}")
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
description: 'Read counts at each stage of the alignment → deduplication → RF-count pipeline per sample.'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-count-progression'
  title: 'nf-core/rnastructurome Count Progression'
  show_table_by_default: true
headers:
  mapped_reads_pre_dedup:
    title: 'Mapped (pre-dedup)'
    description: 'Reads mapped to the reference before deduplication'
    scale: 'Blues'
    format: '{:,.0f}'
  mapped_reads_post_dedup:
    title: 'Mapped (post-dedup)'
    description: 'Reads retained after UMI/positional deduplication'
    scale: 'Blues'
    format: '{:,.0f}'
  pct_removed_by_dedup:
    title: 'Removed by Dedup'
    description: 'Percentage of mapped reads removed as duplicates'
    scale: 'Oranges'
    format: '{:,.1f}'
    suffix: '%'
  rfcount_covered_transcripts:
    title: 'RF-count: Covered Transcripts'
    description: 'Number of transcripts with sufficient coverage in RF-count'
    scale: 'Greens'
    format: '{:,.0f}'
data:
${dataBlock}
"""
}

def parseRfnormLog(logFile) {
    def covered   = 0L
    def discarded = 0L
    logFile.readLines().each { line ->
        def covM = (line =~ /\[\*\]\s+Covered transcripts:\s+(\d+)/)
        if (covM.find()) covered = covM.group(1) as long
        def disM = (line =~ /\[\*\]\s+Discarded transcripts:\s+(\d+)\s+total/)
        if (disM.find()) discarded = disM.group(1) as long
    }
    [covered: covered, discarded: discarded]
}

def parseRffoldLog(logFile) {
    def folded    = 0L
    def discarded = 0L
    logFile.readLines().each { line ->
        def foldM = (line =~ /\[\*\]\s+Folded transcripts:\s+(\d+)/)
        if (foldM.find()) folded = foldM.group(1) as long
        def disM  = (line =~ /\[\*\]\s+Discarded transcripts:\s+(\d+)\s+total/)
        if (disM.find()) discarded = disM.group(1) as long
    }
    [folded: folded, discarded: discarded]
}

def buildSimpleMultiqcTable(rows, id, sectionName, description, headers) {
    def rowEntries
    if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else if (rows instanceof List && rows.every { row -> row instanceof List && row.size() == 2 && row[1] instanceof Map }) {
        rowEntries = rows
    } else {
        rowEntries = []
    }
    def orderedRows = rowEntries.sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sampleId = row[0]
        def metrics  = row[1]
        def metricLines = metrics.collect { key, value -> "    ${key}: ${value}" }.join('\n')
        "  ${sampleId}:\n${metricLines}"
    }.join('\n')
    def headerBlock = headers.collect { col, cfg ->
        def lines = ["  ${col}:"]
        cfg.each { k, v -> lines << "    ${k}: '${v}'" }
        lines.join('\n')
    }.join('\n')

    """id: '${id}'
section_name: '${sectionName}'
description: '${description}'
plot_type: 'table'
pconfig:
  id: '${id}'
  title: '${sectionName}'
  show_table_by_default: true
headers:
${headerBlock}
data:
${dataBlock}
"""
}

def rfnormStatsMultiqc(rows) {
    buildSimpleMultiqcTable(
        rows,
        'nf-core-rnastructurome-rfnorm-stats',
        'nf-core/rnastructurome RF-norm Statistics',
        'Transcript coverage and discard counts from rf-norm (per normalisation group).',
        [
            covered  : [title: 'Covered Transcripts',  description: 'Transcripts with sufficient coverage for normalisation', scale: 'Greens', format: '{:,.0f}'],
            discarded: [title: 'Discarded Transcripts', description: 'Transcripts discarded by rf-norm (insufficient coverage, mismatches, or absent in control)', scale: 'Reds', format: '{:,.0f}']
        ]
    )
}

def rffoldStatsMultiqc(rows) {
    buildSimpleMultiqcTable(
        rows,
        'nf-core-rnastructurome-rffold-stats',
        'nf-core/rnastructurome RF-fold Statistics',
        'Folded and discarded transcript counts from rf-fold (per fold group).',
        [
            folded   : [title: 'Folded Transcripts',   description: 'Transcripts successfully folded by rf-fold', scale: 'Purples', format: '{:,.0f}'],
            discarded: [title: 'Discarded Transcripts', description: 'Transcripts discarded by rf-fold (XML parse errors or folding failures)', scale: 'Reds', format: '{:,.0f}']
        ]
    )
}

def cutadaptAdaptersMultiqc(rows) {
    def rowEntries
    if (rows instanceof Map) {
        rowEntries = rows.entrySet().collect { entry -> [entry.key, entry.value] }
    } else if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        rowEntries = [rows]
    } else if (rows instanceof List) {
        rowEntries = rows.collectMany { row ->
            if (row instanceof Map) {
                return row.entrySet().collect { entry -> [entry.key, entry.value] }
            }
            if (row instanceof List && row.size() == 2 && row[1] instanceof Map) {
                return [row]
            }
            return []
        }
    } else {
        rowEntries = []
    }
    def orderedRows = rowEntries
        .findResults { row ->
            if (row instanceof List && row.size() == 2 && row[1] instanceof Map) {
                return [row[0].toString(), row[1]]
            }
            null
        }
        .sort { a, b -> a[0] <=> b[0] }
    def dataBlock = orderedRows.collect { row ->
        def sample_id = row[0]
        def metrics = row[1]
        def metricLines = metrics.collect { key, value ->
            def rendered = value.toString().replace("'", "''")
            "    ${key}: '${rendered}'"
        }.join('\n')
        "  ${sample_id}:\n${metricLines}"
    }.join('\n')

    """id: 'nf-core-rnastructurome-cutadapt-adapters'
section_name: 'Cutadapt: Adapter Sequences Used'
description: 'Adapter sequences used for trimming in each sample.'
parent_id: 'cutadapt'
parent_name: 'Cutadapt'
plot_type: 'table'
pconfig:
  id: 'nf-core-rnastructurome-cutadapt-adapters'
  title: 'Cutadapt: Adapter Sequences Used'
  show_table_by_default: true
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
