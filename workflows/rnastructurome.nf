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
include { STAR_GENOMEGENERATE                     } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN as STAR_ALIGN_RTSTOP         } from '../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_MAP            } from '../modules/nf-core/star/align/main'
include { BOWTIE_BUILD          } from '../modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN          } from '../modules/nf-core/bowtie/align/main'
include { BOWTIE2_BUILD         } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN         } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT                              } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FIXMATE      } from '../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SORT  } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MARKDUP      } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_STATS        } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT     } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_PRE } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS     } from '../modules/nf-core/samtools/idxstats/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { RNAFRAMEWORK_RFCOUNT         } from '../modules/local/rnaframework/count/main'
include { RNAFRAMEWORK_RFCOUNT_GENOME      } from '../modules/local/rnaframework/count_genome/main'
include { RNAFRAMEWORK_RFRCTOOLS_EXTRACT   } from '../modules/local/rnaframework/rctools/extract/main'
include { BEDOPS_GTF2BED               } from '../modules/nf-core/bedops/gtf2bed/main'
include { RSEQC_INFEREXPERIMENT        } from '../modules/nf-core/rseqc/inferexperiment/main'
include { RNAFRAMEWORK_RFNORM      } from '../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFJACKKNIFE } from '../modules/local/rnaframework/jackknife/main'
include { RNAFRAMEWORK_RFFOLD      } from '../modules/local/rnaframework/fold/main'
include { ENSEMBL_TRANSCRIPTOME } from '../modules/local/ensembl/transcriptome/main'
include { ENSEMBL_GENOME        } from '../modules/local/ensembl/genome/main'
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
include { WIG_TO_GENOME as WIG_TO_GENOME_REACTIVITY } from '../modules/local/wig_to_genome/main'
include { WIG_TO_GENOME as WIG_TO_GENOME_SHANNON    } from '../modules/local/wig_to_genome/main'
include { RNAFRAMEWORK_TORDAT    } from '../modules/local/tordat/main'
include { R2DT                   } from '../modules/local/r2dt/main'
include { VIENNARNA              } from '../modules/local/viennarna/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_REACTIVITY            } from '../modules/nf-core/ucsc/wigtobigwig/main'
include { UCSC_WIGTOBIGWIG as UCSC_WIGTOBIGWIG_SHANNON               } from '../modules/nf-core/ucsc/wigtobigwig/main'

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
            // User-supplied local FASTA: genome_fasta for STAR route, transcriptome_fasta for
            // --transcriptome route.  The legacy --fasta flag maps to the appropriate route.
            def local_fasta = pipeline_config.transcriptome
                ? (pipeline_config.transcriptome_fasta ?: pipeline_config.fasta)
                : (pipeline_config.genome_fasta        ?: pipeline_config.fasta)
            if (local_fasta) {
                return [ reference_key, "path::${local_fasta.toString()}", original_organism ]
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

    // Ensembl downloads are mutually exclusive: genome mode downloads the soft-masked
    // genome FASTA (for STAR); transcriptome mode downloads the cDNA FASTA (for Bowtie).
    // ch_ensembl_not_found triggers NCBI fallback for species absent from Ensembl in both modes.
    // ch_ensembl_fasta_source_url carries the download URL(s) for provenance reporting.
    def ch_ensembl_not_found        = channel.empty()
    def ch_ensembl_fasta_source_url = channel.empty()
    def ch_reference_genome_fasta_keyed = channel.empty()

    if (pipeline_config.transcriptome) {
        ENSEMBL_TRANSCRIPTOME (
            ch_reference_ensembl_input,
            [
                ensembl_release : pipeline_config.ensembl_release,
                ensembl_base_url: pipeline_config.ensembl_base_url
            ],
            file("${projectDir}/bin/ensembl_transcriptome.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(ENSEMBL_TRANSCRIPTOME.out.versions)
        ch_ensembl_not_found        = ENSEMBL_TRANSCRIPTOME.out.not_found
        ch_ensembl_fasta_source_url = ENSEMBL_TRANSCRIPTOME.out.source_urls
    }

    // MODULE: ENSEMBL_GENOME — download soft-masked genome FASTA for STAR alignment.
    if (!pipeline_config.transcriptome) {
        ENSEMBL_GENOME(
            ch_reference_ensembl_input,
            [
                ensembl_release : pipeline_config.ensembl_release,
                ensembl_base_url: pipeline_config.ensembl_base_url
            ],
            file("${projectDir}/bin/ensembl_genome.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(ENSEMBL_GENOME.out.versions)
        ch_reference_genome_fasta_keyed = ENSEMBL_GENOME.out.fasta
            .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
        ch_ensembl_not_found        = ENSEMBL_GENOME.out.not_found
        ch_ensembl_fasta_source_url = ENSEMBL_GENOME.out.source_url
    }

    // Organisms not found on Ensembl FTP are routed to NCBI_FASTA for automatic accession
    // search.  meta.original_organism carries the raw samplesheet organism string used as
    // the esearch query when no accessions are pre-configured.
    ch_reference_ncbi_from_ensembl = ch_ensembl_not_found
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

    // Keyed GTF channel for .join() in STAR alignment blocks.
    // Format: [ ref_key, meta, gtf ] — allows join with keyed sample channels.
    def ch_gtf_for_star_join = ch_all_reference_gtf.map { meta, gtf ->
        [ meta.id.toString(), meta, gtf ]
    }

    def ch_fasta_sort_script = file("${projectDir}/bin/fasta_sort.py", checkIfExists: true)

    FASTA_SORT_LOCAL (
        ch_reference_local,
        ch_fasta_sort_script
    )
    ch_versions = ch_versions.mix(FASTA_SORT_LOCAL.out.versions)

    def ch_fasta_sort_ensembl_out = channel.empty()
    if (pipeline_config.transcriptome) {
        FASTA_SORT_ENSEMBL (
            ENSEMBL_TRANSCRIPTOME.out.fasta,
            ch_fasta_sort_script
        )
        ch_versions = ch_versions.mix(FASTA_SORT_ENSEMBL.out.versions)
        ch_fasta_sort_ensembl_out = FASTA_SORT_ENSEMBL.out.fasta
    }

    FASTA_SORT_NCBI (
        NCBI_FASTA.out.fasta,
        ch_fasta_sort_script
    )
    ch_versions = ch_versions.mix(FASTA_SORT_NCBI.out.versions)

    ch_reference_fasta_keyed = FASTA_SORT_LOCAL.out.fasta
        .mix(ch_fasta_sort_ensembl_out)
        .mix(FASTA_SORT_NCBI.out.fasta)
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    ch_reference_gtf_keyed   = ch_all_reference_gtf.map { meta, gtf -> [ meta.id.toString(), [meta, gtf] ] }
    ch_reference_fasta_map   = ch_reference_fasta_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
        .first()
    ch_reference_gtf_map     = ch_reference_gtf_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
        .first()

    // Genome FASTA map for STAR index building.
    // For Ensembl species: ENSEMBL_GENOME output (soft-masked genome).
    // For NCBI species (bacteria, viruses): NCBI_FASTA output serves as the genome
    // reference (no introns — genome and transcriptome are equivalent).
    // For local FASTA references on the STAR route: the user-supplied FASTA is
    // the genome (e.g. a viral/mitochondrial genome or a custom assembly).
    def ch_reference_genome_ncbi_keyed = NCBI_FASTA.out.fasta
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    ch_reference_genome_fasta_keyed = ch_reference_genome_fasta_keyed
        .mix(ch_reference_genome_ncbi_keyed)
    if (!pipeline_config.transcriptome) {
        ch_reference_genome_fasta_keyed = ch_reference_genome_fasta_keyed
            .mix(FASTA_SORT_LOCAL.out.fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] })
    }

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
    // INDEX BUILDING — conditional on chosen aligner per principle
    //
    def ch_star_index_map    = channel.value([:])
    def ch_star_index_keyed  = channel.empty()
    def ch_bowtie_index_map  = channel.value([:])
    def ch_bowtie2_index_map = channel.value([:])

    if (!pipeline_config.transcriptome) {
        // Build one STAR index per reference using the genome FASTA + GTF.
        // For Ensembl species this is the soft-masked toplevel genome assembly.
        // For NCBI species (bacteria, viruses) the NCBI FASTA is used as the genome.
        def ch_star_build = ch_reference_genome_fasta_keyed
            .join(ch_reference_gtf_keyed, remainder: true)
            .map { key, fasta_tuple, gtf_tuple ->
                def fasta_meta = fasta_tuple[0]
                def fasta      = fasta_tuple[1]
                def gtf_meta   = gtf_tuple ? gtf_tuple[0] : [id: "${key}_gtf"]
                def gtf        = gtf_tuple ? gtf_tuple[1] : []
                [ [fasta_meta, fasta], [gtf_meta, gtf] ]
            }

        STAR_GENOMEGENERATE(
            ch_star_build.map { entry -> entry[0] },
            ch_star_build.map { entry -> entry[1] }
        )

        ch_star_index_map = STAR_GENOMEGENERATE.out.index
            .map { meta, index -> [ (meta.id.toString()): [meta, index] ] }
            .collect()
            .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
            .first()
        ch_star_index_keyed = STAR_GENOMEGENERATE.out.index
            .map { meta, index -> [ meta.id.toString(), meta, index ] }
    }

    if (pipeline_config.transcriptome) {
        BOWTIE_BUILD(ch_rtstop_reference_fasta)
        ch_bowtie_index_map = BOWTIE_BUILD.out.index
            .map { meta, index -> [ (meta.id.toString()): [meta, index] ] }
            .collect()
            .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
            .first()
    }

    if (pipeline_config.transcriptome) {
        BOWTIE2_BUILD(ch_map_reference_fasta)
        ch_bowtie2_index_map = BOWTIE2_BUILD.out.index
            .map { meta, index -> [ (meta.id.toString()): [meta, index] ] }
            .collect()
            .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
            .first()
    }

    //
    // RT-STOP ALIGNMENT
    //
    def ch_rtstop_aligned_bam      = channel.empty()

    if (!pipeline_config.transcriptome) {
        def ch_rtstop_keyed = ch_rtstop_trimmed_for_align
            .map { meta, reads -> [ resolveReferenceKey(meta, pipeline_config.organism), meta, reads ] }
        def ch_rtstop_star_inputs = ch_rtstop_keyed
            .join(ch_star_index_keyed)
            .join(ch_gtf_for_star_join, remainder: true)
            .map { ref_key, sample_meta, reads, idx_meta, index, gtf_meta, gtf ->
                def has_gtf = gtf_meta != null
                [ [sample_meta, reads], [idx_meta, index], [gtf_meta ?: [id: 'no_gtf'], gtf ?: []], !has_gtf ]
            }
        def ch_rtstop_star_split = ch_rtstop_star_inputs.multiMap { entry ->
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
        def ch_map_keyed = ch_map_trimmed_for_align
            .map { meta, reads -> [ resolveReferenceKey(meta, pipeline_config.organism), meta, reads ] }
        def ch_map_star_inputs = ch_map_keyed
            .join(ch_star_index_keyed)
            .join(ch_gtf_for_star_join, remainder: true)
            .map { ref_key, sample_meta, reads, idx_meta, index, gtf_meta, gtf ->
                def has_gtf = gtf_meta != null
                [ [sample_meta, reads], [idx_meta, index], [gtf_meta ?: [id: 'no_gtf'], gtf ?: []], !has_gtf ]
            }
        def ch_map_star_split = ch_map_star_inputs.multiMap { entry ->
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
    // MODULES: BEDOPS_GTF2BED + RSEQC_INFEREXPERIMENT (STAR route only)
    // Convert each reference GTF to BED12 once, then run infer_experiment on
    // every final BAM to determine library strandedness automatically.
    //
    def ch_strandedness_by_id = channel.empty()

    if (!pipeline_config.transcriptome) {
        BEDOPS_GTF2BED(ch_all_reference_gtf)

        def ch_reference_bed_map = BEDOPS_GTF2BED.out.bed
            .map { meta, bed -> [ (meta.id.toString()): [meta, bed] ] }
            .collect()
            .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }

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

    //
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
        def ch_genome_fasta_map = ch_reference_genome_fasta_keyed
            .map { key, value -> [ (key): value ] }
            .collect()
            .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }

        def ch_rfcount_genome_inputs = ch_markdup_bam_bai
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
        ch_rfcount_summary = RNAFRAMEWORK_RFCOUNT_GENOME.out.summary
        ch_rfcount_plots   = RNAFRAMEWORK_RFCOUNT_GENOME.out.plots
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFCOUNT_GENOME.out.versions)

        // rf-rctools extract: genome RC → transcript-level RC using reference GTF.
        // The module generates per-file .rci indexes itself (rf-rctools index) and calls
        // rf-rctools extract with the BASENAME so strand-aware extraction is activated.
        def ch_rctools_inputs = RNAFRAMEWORK_RFCOUNT_GENOME.out.rc
            .combine(ch_reference_gtf_map)
            .map { combined ->
                def meta    = combined[0]
                def rc      = combined[1]
                def gtf_map = combined[2]
                def ref_key = resolveReferenceKey(meta, pipeline_config.organism)
                def gtf_t   = gtf_map[ref_key]
                if (!gtf_t) error("No GTF resolved for reference '${ref_key}' for rf-rctools extract.")
                [ [meta, rc, []], gtf_t ]
            }
        def ch_rct_split = ch_rctools_inputs.multiMap { entry ->
            rc:  entry[0]
            gtf: entry[1]
        }
        RNAFRAMEWORK_RFRCTOOLS_EXTRACT(ch_rct_split.rc, ch_rct_split.gtf)
        ch_rfcount_rc  = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rc
        ch_rfcount_rci = RNAFRAMEWORK_RFRCTOOLS_EXTRACT.out.rci
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

    def ch_pre_dedup_mapped_reads = SAMTOOLS_FLAGSTAT_PRE.out.flagstat
        .map { meta, flagstat -> [ meta.id.toString(), parseFlagstatMappedReads(flagstat) ] }

    def ch_post_dedup_mapped_reads = SAMTOOLS_FLAGSTAT.out.flagstat
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

    //
    // MODULE: merge_bp — merge per-transcript .bp files into a single file per fold group for genome browser visualisation
    //
    MERGE_BP (
        RNAFRAMEWORK_DOTPLOT2BP.out.bp,
        file("${projectDir}/bin/merge_bp.py", checkIfExists: true)
    )

    //
    // MODULE: R2DT — template-based 2D structure diagrams with reactivity overlay
    //
    if (params.r2dt) {
        def ch_r2dt_xml = RNAFRAMEWORK_RFNORM.out.xml
            .map { meta, xml ->
                def fold_group = meta.cell_line?.toString()
                if (!fold_group) error("Missing cell_line for '${meta.id}' — required for R2DT grouping.")
                [ fold_group, xml instanceof List ? xml : [xml] ]
            }
            .groupTuple()
            .map { fold_group, xml_lists -> [ fold_group, xml_lists.flatten() ] }

        def ch_r2dt_input = RNAFRAMEWORK_RFFOLD.out.structures
            .map { meta, dir -> [ meta.id.toString(), meta, dir ] }
            .join(ch_r2dt_xml)
            .map { _fg, meta, dir, xmls -> [ meta, dir, xmls ] }
            .combine(ch_reference_fasta_map)
            .flatMap { combined ->
                def meta      = combined[0]
                def fold_dir  = combined[1]
                def xmls      = combined[2]
                def fasta_map = combined[3]
                def ref_key   = resolveReferenceKey(meta, pipeline_config.organism)
                def fasta_t   = fasta_map[ref_key]
                if (!fasta_t) {
                    log.warn("Skipping R2DT for '${meta.id}': no FASTA for '${ref_key}'")
                    return []
                }
                return [ [ meta, fold_dir, xmls, fasta_t[1] ] ]
            }

        R2DT(
            ch_r2dt_input,
            file("${projectDir}/bin/r2dt_colour_svg.py",         checkIfExists: true),
            file("${projectDir}/bin/r2dt_extract_sequences.py",  checkIfExists: true)
        )
        ch_versions = ch_versions.mix(R2DT.out.versions.first())

        def ch_rnaplot_input = RNAFRAMEWORK_RFFOLD.out.structures
            .map { meta, dir -> [ meta.id.toString(), meta, dir ] }
            .join(ch_fold_input.map { meta, xmls -> [ meta.id.toString(), xmls ] })
            .join(R2DT.out.drawn_ids.map { meta, f -> [ meta.id.toString(), f ] })
            .map { _id, meta, dir, xmls, drawn -> [ meta, dir, xmls, drawn ] }

        VIENNARNA(
            ch_rnaplot_input,
            file("${projectDir}/bin/viennarna_extract_xml.py",   checkIfExists: true),
            file("${projectDir}/bin/viennarna_colour_svg.py",    checkIfExists: true)
        )
        ch_versions = ch_versions.mix(VIENNARNA.out.versions.first())
    } else {
        // R2DT is container-only; when running without containers draw all structures with ViennaRNA
        def ch_rnaplot_input = RNAFRAMEWORK_RFFOLD.out.structures
            .map { meta, dir -> [ meta.id.toString(), meta, dir ] }
            .join(ch_fold_input.map { meta, xmls -> [ meta.id.toString(), xmls ] })
            .map { _id, meta, dir, xmls -> [ meta, dir, xmls, file('/dev/null') ] }

        VIENNARNA(
            ch_rnaplot_input,
            file("${projectDir}/bin/viennarna_extract_xml.py",   checkIfExists: true),
            file("${projectDir}/bin/viennarna_colour_svg.py",    checkIfExists: true)
        )
        ch_versions = ch_versions.mix(VIENNARNA.out.versions.first())
    }

    //
    // MODULE: rf-wiggle — convert rf-norm XML reactivities to WIG + chrom.sizes
    //
    RNAFRAMEWORK_RFWIGGLE (
        RNAFRAMEWORK_RFNORM.out.xml
    )

    MERGE_WIG (
        RNAFRAMEWORK_RFWIGGLE.out.wig
    )

    //
    // Group per-replicate merged WIGs by cell_line.
    // Single replicate: bypass AVERAGE_WIG, keep meta.id = "HEK293T_1" → HEK293T_1_reactivity.bw
    // Multiple replicates: run AVERAGE_WIG, set meta.id = cell_line   → HEK293T_reactivity.bw
    //
    def ch_reactivity_grouped = MERGE_WIG.out.merged_wig
        .map { meta, wig -> [ meta.cell_line.toString(), meta, wig ] }
        .groupTuple(by: 0)
        .map { cell_line, metas, wigs ->
            def base_meta = wigs.size() == 1
                ? metas[0]
                : metas[0] + [ id: cell_line ]
            [ base_meta, wigs.flatten() ]
        }

    def ch_reactivity_branches = ch_reactivity_grouped.branch { _meta, wigs ->
        single: wigs.size() == 1
        multi:  true
    }

    AVERAGE_WIG (
        ch_reactivity_branches.multi,
        file("${projectDir}/bin/average_wig.py", checkIfExists: true)
    )

    //
    // MODULE: wig_to_genome + wigToBigWig — convert merged/averaged transcript WIG
    // to genomic BigWig for IGV
    //
    def ch_reactivity_for_bigwig = ch_reactivity_branches.single
        .map { meta, wigs -> [ meta, wigs[0] ] }
        .mix(AVERAGE_WIG.out.merged_wig)
    def ch_reactivity_genomic_wig_input = ch_reactivity_for_bigwig
        .combine(ch_reference_gtf_map)
        .flatMap { combined ->
            def meta = combined[0]
            def wig = combined[1]
            def gtf_map = combined[2]
            def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
            def gtf_tuple = gtf_map[reference_key]
            if (!gtf_tuple) {
                log.warn("Skipping genomic reactivity BigWig for '${reference_key}': no GTF available.")
                return []
            }
            return [ [ meta, wig, gtf_tuple[1] ] ]
        }

    WIG_TO_GENOME_REACTIVITY (
        ch_reactivity_genomic_wig_input,
        file("${projectDir}/bin/remap_wig_to_genome.py", checkIfExists: true)
    )

    UCSC_WIGTOBIGWIG_REACTIVITY (
        WIG_TO_GENOME_REACTIVITY.out.wig,
        WIG_TO_GENOME_REACTIVITY.out.chrom_sizes.map { _meta, sizes -> sizes }
    )

    //
    // MODULE: merge_shannon_wig + wigToBigWig — merge per-transcript Shannon entropy WIG files
    // and convert to BigWig for genome browser visualisation
    //
    if (params.rffold_shannon_entropy) {
        MERGE_SHANNON_WIG (
            RNAFRAMEWORK_RFFOLD.out.shannon_wig
        )

        def ch_shannon_genomic_wig_input = MERGE_SHANNON_WIG.out.merged_wig
            .combine(ch_reference_gtf_map)
            .flatMap { combined ->
                def meta = combined[0]
                def wig = combined[1]
                def gtf_map = combined[2]
                def reference_key = resolveReferenceKey(meta, pipeline_config.organism)
                def gtf_tuple = gtf_map[reference_key]
                if (!gtf_tuple) {
                    log.warn("Skipping genomic Shannon BigWig for '${reference_key}': no GTF available.")
                    return []
                }
                return [ [ meta, wig, gtf_tuple[1] ] ]
            }

        WIG_TO_GENOME_SHANNON (
            ch_shannon_genomic_wig_input,
            file("${projectDir}/bin/remap_wig_to_genome.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(WIG_TO_GENOME_SHANNON.out.versions.first())

        UCSC_WIGTOBIGWIG_SHANNON (
            WIG_TO_GENOME_SHANNON.out.wig,
            WIG_TO_GENOME_SHANNON.out.chrom_sizes.map { _meta, sizes -> sizes }
        )
    }

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
        ch_ref_fasta_names = FASTA_SORT_LOCAL.out.fasta
            .map { meta, _fasta -> [ meta.id.toString(), local_fasta_name ] }
        ch_ref_gtf_names = FASTA_SORT_LOCAL.out.fasta
            .map { meta, _fasta -> [ meta.id.toString(), local_gtf_name ] }
    } else {
        ch_ref_fasta_names = ch_ensembl_fasta_source_url
            .map { meta, urls_file ->
                def names = urls_file.readLines().findAll { line -> line.trim() }
                    .collect { line -> line.tokenize('/').last() }.join(' + ')
                [ meta.id.toString(), names ]
            }
            .mix(NCBI_FASTA.out.source_accessions
                .map { meta, acc_file ->
                    def accs = acc_file.readLines()
                        .findAll { line -> line.trim() && !line.startsWith('stub:') }
                        .collect { line -> line.tokenize('/').last() }.join(', ')
                    [ meta.id.toString(), accs ?: meta.id.toString() ]
                })
        ch_ref_gtf_names = ENSEMBL_GTF.out.source_urls
            .map { meta, urls_file ->
                def name = urls_file.readLines().find { line -> line.trim() }?.tokenize('/')?.last() ?: ''
                [ meta.id.toString(), name ]
            }
            .mix(ch_reference_gtf_local
                .map { meta, gtf_file -> [ meta.id.toString(), gtf_file.name ] })
            .mix(NCBI_FASTA.out.source_accessions
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
    ch_versions = ch_versions.mix(MERGE_BP.out.versions.first())
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFWIGGLE.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_WIG.out.versions.first())
    ch_versions = ch_versions.mix(AVERAGE_WIG.out.versions.first())
    ch_versions = ch_versions.mix(WIG_TO_GENOME_REACTIVITY.out.versions.first())
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

def parseInferExperiment(txtFile) {
    def forward = 0.0
    def reverse = 0.0
    txtFile.readLines().each { line ->
        def m = line =~ /Fraction of reads explained by "(?:1\+\+,1--,2\+-,2-\+|\+\+,--)": (.+)/
        if (m) forward = m[0][1].trim() as double
        m = line =~ /Fraction of reads explained by "(?:1\+-,1-\+,2\+\+,2--|\\+-,-\+)": (.+)/
        if (m) reverse = m[0][1].trim() as double
    }
    if (forward > 0.7) return 'first'
    if (reverse > 0.7) return 'second'
    return 'unstranded'
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

    if (pipeline_config.transcriptome) {
        moduleOptions['aligner'] = 'transcriptome'
        if (!principles || principles.contains('rt-stop')) {
            moduleOptions['rtstop_aligner'] = 'bowtie'
            moduleOptions['rtstop_aligner_args'] = renderBowtie1Args(pipeline_config)
        }
        if (principles.contains('map')) {
            moduleOptions['map_aligner'] = 'bowtie2'
            moduleOptions['map_aligner_args'] = renderBowtie2Args(pipeline_config)
        }
    } else {
        moduleOptions['aligner'] = 'star'
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
    if (!(requestedMethod in [2, 3, 4])) {
        error("Unsupported rf-norm normalization method '${pipeline_config.rfnorm_norm_method}'. Expected one of: 2, 3, 4.")
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
        3: 'Box-plot normalisation',
        4: 'Mitchell normalisation'
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
        // Single row, collect() flattened [id, map] to a bare 2-element list
        rowEntries = [rows]
    } else if (rows instanceof List && rows.every { row -> row instanceof List && row.size() == 2 && row[1] instanceof Map }) {
        // Multiple rows as nested [[id, map], ...] (collect(flat:false) behaviour)
        rowEntries = rows
    } else if (rows instanceof List && rows.size() % 2 == 0 && rows.collate(2).every { pair -> pair.size() == 2 && pair[1] instanceof Map }) {
        // Multiple rows flattened by Nextflow collect() into [id1, map1, id2, map2, ...]
        rowEntries = rows.collate(2)
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
        'Transcript coverage statistics from RF-count and RF-norm (per normalisation group).',
        [
            rfcount_covered: [title: 'RF-count Covered', description: 'Transcripts covered by RF-count (input to RF-norm)', scale: 'Blues',  format: '{:,.0f}'],
            covered        : [title: 'RF-norm Covered',  description: 'Transcripts passing RF-norm normalisation (sufficient coverage)', scale: 'Greens', format: '{:,.0f}']
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
