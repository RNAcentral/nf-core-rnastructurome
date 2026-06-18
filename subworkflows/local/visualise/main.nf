//
// VISUALISE_STRUCTURES — render 2D structure diagrams from rf-fold structures.
//
// When R2DT is enabled (container-only), draw template-based diagrams with a
// reactivity overlay and fall back to ViennaRNA for structures R2DT could not
// draw. Otherwise draw everything with ViennaRNA. Produces only published SVGs
// and software versions; nothing here is consumed downstream.
//

include { R2DT                } from '../../../modules/local/r2dt/main'
include { VIENNARNA           } from '../../../modules/local/viennarna/main'
include { resolveReferenceKey } from '../../../workflows/rnastructurome_functions.nf'

workflow VISUALISE_STRUCTURES {

    take:
    ch_rfnorm_xml          // channel: [ val(meta), path(xml) ]
    ch_fold_structures     // channel: [ val(meta), path(dir) ]
    ch_fold_input          // channel: [ val(meta), path(xmls) ]
    ch_reference_fasta_map // value:   map of reference_key -> [meta, fasta]
    pipeline_config        // map

    main:
    ch_versions = channel.empty()

    if (params.r2dt) {
        def ch_r2dt_xml = ch_rfnorm_xml
            .map { meta, xml ->
                def fold_group = meta.cell_line?.toString()
                if (!fold_group) error("Missing cell_line for '${meta.id}' — required for R2DT grouping.")
                [ fold_group, xml instanceof List ? xml : [xml] ]
            }
            .groupTuple()
            .map { fold_group, xml_lists -> [ fold_group, xml_lists.flatten() ] }

        def ch_r2dt_input = ch_fold_structures
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

        def ch_rnaplot_input = ch_fold_structures
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
        def ch_rnaplot_input = ch_fold_structures
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

    emit:
    versions = ch_versions
}
