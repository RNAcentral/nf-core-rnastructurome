// PREPARE_REFERENCES — resolve each sample's reference to a local path, Ensembl species, or NCBI
// accessions (with Ensembl-404 fallback), fetch/sort FASTA+GTF into lookup maps, then build the aligner
// index (STAR or Bowtie/Bowtie2). Source-provenance channels are emitted for RDAT naming downstream.

include { ENSEMBL_TRANSCRIPTOME } from '../../../modules/local/ensembl/transcriptome/main'
include { ENSEMBL_GENOME        } from '../../../modules/local/ensembl/genome/main'
include { ENSEMBL_GTF           } from '../../../modules/local/ensembl/gtf/main'
include { NCBI_FASTA            } from '../../../modules/local/ncbi/fasta/main'
include { NCBI_GTF              } from '../../../modules/local/ncbi/gtf/main'
include { FASTA_SORT as FASTA_SORT_LOCAL   } from '../../../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_ENSEMBL } from '../../../modules/local/fasta/sort/main'
include { FASTA_SORT as FASTA_SORT_NCBI    } from '../../../modules/local/fasta/sort/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'
include { BOWTIE_BUILD          } from '../../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD         } from '../../../modules/nf-core/bowtie2/build/main'

include { resolveReferenceResolution } from '../../../workflows/rnastructurome_functions.nf'
include { uniqueReferenceResolution  } from '../../../workflows/rnastructurome_functions.nf'
include { collectToMap               } from '../../../workflows/rnastructurome_functions.nf'
include { resolveReferenceKey        } from '../../../workflows/rnastructurome_functions.nf'

workflow PREPARE_REFERENCES {

    take:
    ch_samplesheet_for_branching // channel: [ val(meta), [ reads ] ] (post-cat)
    pipeline_config              // map

    main:
    ch_versions = channel.empty()

    // Branch by probing principle (for transcriptome-route reference-FASTA selection)
    def principle_branches = ch_samplesheet_for_branching.branch { meta, _reads ->
        rtstop: (meta.principle ?: '').toLowerCase() == 'rt-stop'
        map:    (meta.principle ?: '').toLowerCase() == 'map'
    }

    def ch_reference_requests = uniqueReferenceResolution(
        ch_samplesheet_for_branching.map { meta, _reads -> resolveReferenceResolution(meta, pipeline_config, 'fasta') },
        'transcript reference'
    )

    def ch_reference_local = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('path::') }
        .map { reference_key, resolution, _original_organism ->
            def fasta_path = resolution - 'path::'
            [ [ id: reference_key, organism: reference_key ], file(fasta_path, checkIfExists: true) ]
        }

    def ch_reference_ensembl_input = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ensembl::') }
        .map { reference_key, resolution, original_organism ->
            def ensembl_species = resolution - 'ensembl::'
            [ [ id: reference_key, organism: reference_key, ensembl_species: ensembl_species,
                original_organism: original_organism ], ensembl_species ]
        }

    // Ensembl downloads are mutually exclusive: genome mode gets the soft-masked genome FASTA (STAR),
    // transcriptome mode gets the cDNA FASTA (Bowtie). ch_ensembl_not_found triggers NCBI fallback.
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
    // ch_reference_genome_fasta_keyed is populated after sorting (see FASTA_SORT_ENSEMBL block below).
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
        ch_ensembl_not_found        = ENSEMBL_GENOME.out.not_found
        ch_ensembl_fasta_source_url = ENSEMBL_GENOME.out.source_url
    }

    // Organisms not found on Ensembl FTP are routed to NCBI_FASTA for automatic accession search;
    // meta.original_organism carries the esearch query when no accessions are pre-configured.
    def ch_reference_ncbi_from_ensembl = ch_ensembl_not_found
        .map { meta, _not_found_file -> [ meta, "" ] }

    def ch_reference_ncbi_explicit = ch_reference_requests
        .filter { _reference_key, resolution, _org -> resolution.startsWith('ncbi::') }
        .map { reference_key, resolution, original_organism ->
            def accessions = resolution - 'ncbi::'
            [ [ id: reference_key, organism: reference_key, original_organism: original_organism ], accessions ]
        }

    def ch_reference_ncbi_input = ch_reference_ncbi_explicit.mix(ch_reference_ncbi_from_ensembl)

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

    // GTF resolution is skipped entirely when stop_after_jackknife + transcriptome, since every
    // GTF-consuming step is gated behind one of those, avoiding wasted Ensembl lookups.
    def ch_all_reference_gtf        = channel.empty()
    def ch_reference_gtf_local      = channel.empty()
    def ch_ensembl_gtf_source_urls  = channel.empty()

    if (!(pipeline_config.stop_after_jackknife && pipeline_config.transcriptome)) {
        def ch_reference_gtf_requests = uniqueReferenceResolution(
            ch_samplesheet_for_branching.map { meta, _reads -> resolveReferenceResolution(meta, pipeline_config, 'gtf') },
            'GTF reference'
        )

        // Local GTF (already in GTF format) — used directly
        ch_reference_gtf_local = ch_reference_gtf_requests
            .filter { _k, resolution, _o -> resolution.startsWith('path::') }
            .map { reference_key, resolution, _original_organism ->
                [ [ id: reference_key, organism: reference_key ], file(resolution - 'path::', checkIfExists: true) ]
            }

        def ch_reference_gtf_ensembl_input = ch_reference_gtf_requests
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
        ch_ensembl_gtf_source_urls = ENSEMBL_GTF.out.source_urls

        // NCBI references get their annotation from NCBI_GTF's synthetic GTF (already run above);
        // ENSEMBL_GTF.out.not_found is silently dropped since those organisms have a GTF via NCBI_GTF.
        ch_all_reference_gtf = ch_reference_gtf_local
            .mix(ENSEMBL_GTF.out.gtf)
            .mix(NCBI_GTF.out.gtf)
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
    } else {
        // Genome route: sort the Ensembl genome FASTA so it is published to reference/
        // and so STAR_GENOMEGENERATE receives a chromosome-sorted FASTA.
        FASTA_SORT_ENSEMBL (
            ENSEMBL_GENOME.out.fasta,
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

    def ch_reference_fasta_keyed = FASTA_SORT_LOCAL.out.fasta
        .mix(ch_fasta_sort_ensembl_out)
        .mix(FASTA_SORT_NCBI.out.fasta)
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    def ch_reference_gtf_keyed   = ch_all_reference_gtf.map { meta, gtf -> [ meta.id.toString(), [meta, gtf] ] }
    def ch_reference_fasta_map   = collectToMap(ch_reference_fasta_keyed)
    def ch_reference_gtf_map     = collectToMap(ch_reference_gtf_keyed)

    // Genome FASTA map for STAR index building: Ensembl species use ENSEMBL_GENOME's soft-masked genome;
    // NCBI species (no introns) use the sorted NCBI FASTA (must be sorted/uncompressed — STAR/samtools
    // faidx reject the raw plain-gzip download); local FASTA references use the user-supplied FASTA as-is.
    def ch_reference_genome_ncbi_keyed = FASTA_SORT_NCBI.out.fasta
        .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
    if (!pipeline_config.transcriptome) {
        // Genome route: use the sorted Ensembl genome FASTA for STAR index building.
        ch_reference_genome_fasta_keyed = FASTA_SORT_ENSEMBL.out.fasta
            .map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] }
            .mix(ch_reference_genome_ncbi_keyed)
            .mix(FASTA_SORT_LOCAL.out.fasta.map { meta, fasta -> [ meta.id.toString(), [meta, fasta] ] })
    } else {
        ch_reference_genome_fasta_keyed = ch_reference_genome_fasta_keyed
            .mix(ch_reference_genome_ncbi_keyed)
    }

    def ch_rtstop_reference_fasta = principle_branches.rtstop
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

    def ch_map_reference_fasta = principle_branches.map
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

    // INDEX BUILDING — conditional on chosen aligner per principle
    def ch_star_index        = channel.empty()
    def ch_bowtie_index_map  = channel.value([:])
    def ch_bowtie2_index_map = channel.value([:])

    if (!pipeline_config.transcriptome) {
        // Build one STAR index per reference using the genome FASTA + GTF (Ensembl soft-masked
        // toplevel assembly, or the NCBI FASTA for bacteria/viruses).
        def ch_star_build = ch_reference_genome_fasta_keyed
            .join(ch_reference_gtf_keyed, remainder: true)
            .map { key, fasta_tuple, gtf_tuple ->
                def fasta_meta = fasta_tuple[0]
                def fasta      = fasta_tuple[1]
                def gtf_meta   = gtf_tuple ? gtf_tuple[0] : [id: "${key}_gtf"]
                def gtf        = gtf_tuple ? gtf_tuple[1] : []
                [ [fasta_meta, fasta], [gtf_meta, gtf] ]
            }

        def ch_star_build_split = ch_star_build.multiMap { fasta_entry, gtf_entry ->
            fasta: fasta_entry
            gtf:   gtf_entry
        }
        STAR_GENOMEGENERATE(ch_star_build_split.fasta, ch_star_build_split.gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.index
    }

    if (pipeline_config.transcriptome) {
        BOWTIE_BUILD(ch_rtstop_reference_fasta)
        ch_bowtie_index_map = collectToMap(
            BOWTIE_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
        )
    }

    if (pipeline_config.transcriptome) {
        BOWTIE2_BUILD(ch_map_reference_fasta)
        ch_bowtie2_index_map = collectToMap(
            BOWTIE2_BUILD.out.index.map { meta, index -> [ meta.id.toString(), [meta, index] ] }
        )
    }

    emit:
    fasta_map               = ch_reference_fasta_map          // value: map ref_key -> [meta, fasta]
    gtf_map                 = ch_reference_gtf_map            // value: map ref_key -> [meta, gtf]
    all_gtf                 = ch_all_reference_gtf            // channel: [ val(meta), path(gtf) ]
    genome_fasta_keyed      = ch_reference_genome_fasta_keyed // channel: [ key, [meta, fasta] ]
    star_index              = ch_star_index                  // channel: [ val(meta), path(index) ] (genome route)
    bowtie_index_map        = ch_bowtie_index_map            // value: map ref_key -> [meta, index]
    bowtie2_index_map       = ch_bowtie2_index_map           // value: map ref_key -> [meta, index]
    ensembl_fasta_source_url = ch_ensembl_fasta_source_url   // channel: [ val(meta), path(urls) ]
    gtf_local               = ch_reference_gtf_local         // channel: [ val(meta), path(gtf) ]
    local_fasta_sorted      = FASTA_SORT_LOCAL.out.fasta     // channel: [ val(meta), path(fasta) ] (--fasta route naming)
    ncbi_source_accessions  = NCBI_FASTA.out.source_accessions // channel: [ val(meta), path(acc) ]
    ensembl_gtf_source_urls = ch_ensembl_gtf_source_urls     // channel: [ val(meta), path(urls) ] — empty when stop_after_jackknife + transcriptome
    versions                = ch_versions
}
