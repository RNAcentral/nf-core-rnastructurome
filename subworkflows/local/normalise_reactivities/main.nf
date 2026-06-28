//
// NORMALISE_REACTIVITIES — group RC files by sample_group/replicate/condition, pair treated
// with untreated controls (exact match, or fuzzy base-token fallback when enabled), then
// run rf-norm to produce per-base reactivity XML files.
//

include { RNAFRAMEWORK_RFNORM         } from '../../../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFNORMFACTOR   } from '../../../modules/local/rnaframework/normfactor/main'
include { RNAFRAMEWORK_RFRCTOOLS_SPLIT } from '../../../modules/local/rnaframework/rctools/split/main'
include {
    sampleGroupBaseToken
    resolveRfNormNormMethod
    resolveReferenceKey
} from '../../../workflows/rnastructurome_functions.nf'

workflow NORMALISE_REACTIVITIES {

    take:
    ch_rfcount_rc        // channel: [ val(meta), path(rc) ]
    ch_rfcount_rci       // channel: [ val(meta), path(rci) ]
    ch_reference_gtf_map // value:   map of reference_key -> [meta, gtf]
    pipeline_config      // map

    main:
    ch_versions = channel.empty()

    // Attach any available .rci sidecar to the RC file for the same sample.
    def ch_rc_with_rci = ch_rfcount_rc
        .map { meta, rc -> [ meta.id.toString(), meta, rc ] }
        .join(
            ch_rfcount_rci
                .map { meta, rci -> [ meta.id.toString(), rci ] },
            remainder: true
        )
        .map { _sample_id, meta, rc, rci -> [ meta, rc, rci ?: [] ] }

    // Key every sample by its normalisation group (sample_group_replicate) and condition.
    def ch_rc_by_group = ch_rc_with_rci
        .map { meta, rc, rci ->
            if (!meta.sample_group || !meta.replicate) {
                error("Missing sample_group or replicate for sample '${meta.id}'. rf-norm requires both columns to pair samples safely.")
            }
            def group     = "${meta.sample_group}_${meta.replicate}".toString()
            def condition = (meta.condition ?: 'treated').toLowerCase()
            [ group, condition, meta, rc, rci ]
        }

    def ch_treated = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'treated' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }
        .groupTuple()

    def ch_untreated = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'untreated' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }

    def ch_denatured = ch_rc_by_group
        .filter  { _group, condition, _meta, _rc, _rci -> condition == 'denatured' }
        .map     { group, _condition, _meta, rc, _rci -> [ group, rc ] }

    // Stage any available .rci sidecars alongside RC files so rf-norm can auto-discover them.
    def ch_group_rci = ch_rc_by_group
        .filter  { _group, _condition, _meta, _rc, rci -> rci }
        .map     { group, _condition, _meta, _rc, rci -> [ group, rci ] }
        .groupTuple()

    // Enforce rf-norm pairing rules and select the representative meta per group.
    // Defined before the fallback channels so ch_treated_no_untreated can join against it.
    // Note: the denatured-without-untreated check is deferred to ch_norm_input assembly
    // so that fuzzy untreated resolution can supply the missing control first.
    def ch_group_meta = ch_rc_by_group
        .map     { group, condition, meta, _rc, _rci -> [ group, [ condition: condition, meta: meta ] ] }
        .groupTuple()
        .map     { group, entries ->
            def conditions = entries.collect { entry -> entry.condition }.toSet()
            if (!conditions.contains('treated') && conditions.contains('untreated')) {
                def offending = entries
                    .findAll { entry -> entry.condition == 'untreated' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }.sort().join(', ')
                error("rf-norm requires a treated sample for sample_group/replicate group '${group}'. Invalid samples: ${offending}")
            }
            if (conditions.contains('denatured') && !conditions.contains('treated')) {
                def offending = entries
                    .findAll { entry -> entry.condition == 'denatured' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }.sort().join(', ')
                error("rf-norm requires a treated sample for denatured controls in sample_group/replicate group '${group}'. Invalid samples: ${offending}")
            }
            // Always use the treated sample's meta so principle drives scoring-method selection correctly.
            [ group, entries.find { entry -> entry.condition == 'treated' }.meta ]
        }

    // Fuzzy untreated pairing fallback (enabled by default; disable with --fuzzy_untreated_pairing false).
    // When a treated group has no exact sample_group+replicate untreated match, the pipeline
    // falls back to an untreated sample sharing the same sample_group base token (the portion
    // before the first underscore, e.g. "MDA-MB-231" from "MDA-MB-231_MTX") at the same
    // replicate.  The fallback errors if ambiguous; otherwise it warns and proceeds.
    // When disabled, unmatched treated groups proceed without an untreated control
    // (scoring method 2 or 4 instead of 1 or 3).
    def ch_resolved_untreated
    if (pipeline_config.fuzzy_untreated_pairing as Boolean) {
        // Each untreated sample keyed as [base_token, replicate, group, rc] for cross-matching.
        def ch_untreated_for_lookup = ch_rc_by_group
            .filter  { _group, condition, _meta, _rc, _rci -> condition == 'untreated' }
            .map     { group, _condition, meta, rc, _rci ->
                [ sampleGroupBaseToken(meta.sample_group.toString()), meta.replicate.toString(), group, rc ]
            }

        // Treated groups with no direct untreated match — candidates for fuzzy lookup.
        def ch_treated_no_untreated = ch_treated
            .join(ch_group_meta)
            .join(ch_untreated, remainder: true)
            .filter { _group, _treated_rcs, _base_meta, untreated_rc -> !untreated_rc }
            .map    { group, _treated_rcs, base_meta, _untreated_rc ->
                [ sampleGroupBaseToken(base_meta.sample_group.toString()), base_meta.replicate.toString(), group ]
            }

        // Cross-product treated-without-untreated × available untreated, filter on base+rep match,
        // then group by treated group to validate uniqueness before selecting the fallback.
        // When ch_untreated_for_lookup is empty (all-treated datasets), this channel is also
        // empty and the remainder:true join below correctly produces null for untreated.
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
                    error("Ambiguous untreated fallback for '${group}': multiple untreated groups share the same sample_group base and replicate: ${candidateGroups}. Use --fuzzy_untreated_pairing false to disable fuzzy matching.")
                }
                log.warn "No exact untreated match for '${group}' — falling back to '${candidates[0].untreated_group}' (shared sample_group base token at same replicate). Set --fuzzy_untreated_pairing false to require exact matches."
                [ group, candidates[0].rc ]
            }

        // Resolved untreated: direct exact match OR fuzzy fallback.
        ch_resolved_untreated = ch_untreated.mix(ch_fallback_untreated)
    } else {
        // Strict mode: only exact sample_group+replicate matches are used.
        // Groups with no exact untreated match proceed without a negative control.
        ch_resolved_untreated = ch_untreated
    }

    def ch_norm_input = ch_treated
        .join(ch_group_meta)
        .join(ch_resolved_untreated, remainder: true)
        .join(ch_denatured, remainder: true)
        .join(ch_group_rci, remainder: true)
        .map { group, treated_rcs, base_meta, untreated_rc, denatured_rc, rci_files ->
            def hasUntreated  = untreated_rc ? true : false
            def hasDenatured  = denatured_rc ? true : false
            // Deferred from ch_group_meta: checked here so fuzzy pairing can supply the untreated first.
            if (hasDenatured && !hasUntreated) {
                error("rf-norm requires an untreated sample for denatured controls in group '${group}'. No untreated was available (neither an exact sample_group+replicate match nor a fuzzy base-token fallback). Add an untreated sample or set --fuzzy_untreated_pairing false.")
            }
            def principle     = (base_meta.principle ?: '').toLowerCase()
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

    // Cross-experiment normalisation via rf-normfactor. Derives one set of transcriptome-wide
    // normalisation factors per reference (across that reference's treated samples, with their matched
    // untreated/denatured), then feeds the factor file to every group's rf-norm via -nf — putting
    // reactivities on a common scale, unlike per-sample box-plot.
    //
    // Enablement is decided PER REFERENCE:
    //   --rfnorm_use_normfactor true   -> force on for every reference
    //   --rfnorm_use_normfactor false  -> force off (always per-sample box-plot)
    //   --rfnorm_use_normfactor null   -> AUTO (default): on for a reference that has paired untreated
    //                                     controls OR more than one treated sample.
    // References that stay off keep an empty factor slot and normalise each group independently.
    def nfRaw          = pipeline_config.rfnorm_use_normfactor
    def nfForceOn      = (nfRaw != null) && (nfRaw.toString().toLowerCase() in ['true', '1', 'yes'])
    def nfForceOff     = (nfRaw != null) && (nfRaw.toString().toLowerCase() in ['false', '0', 'no'])
    def chunkingActive = pipeline_config.rfnorm_chunk_size && !pipeline_config.transcriptome
    if (nfForceOn && chunkingActive) {
        error("--rfnorm_use_normfactor true is incompatible with --rfnorm_chunk_size: chunking would fragment the cross-experiment factor calculation. Disable one of them.")
    }

    def ch_norm_input_final
    if (nfForceOff || chunkingActive) {
        // rf-normfactor disabled (explicit off, or the chunking path owns normalisation).
        ch_norm_input_final = ch_norm_input
            .map { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files ->
                [ gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, [] ]
            }
    } else {
        // Build per-reference candidates from ch_norm_input, which already pairs each group's treated
        // sample(s) with their RESOLVED untreated/denatured (exact or fuzzy). rf-normfactor pairs
        // -t/-u/-d positionally (treated[i] ↔ untreated[i]), so emit one entry per treated RC carrying
        // its own controls — keeping the pairing intact regardless of groupTuple ordering.
        def ch_nf_pairs = ch_norm_input
            .flatMap { gmeta, treated_rcs, untreated_rc, denatured_rc, _rci ->
                def ref   = resolveReferenceKey(gmeta, pipeline_config.organism)
                def tlist = treated_rcs instanceof List ? treated_rcs : [treated_rcs]
                // A single resolved untreated/denatured control is shared across the group's treated RCs.
                def untreated = (untreated_rc instanceof List ? (untreated_rc ? untreated_rc[0] : null) : untreated_rc) ?: null
                def denatured = (denatured_rc instanceof List ? (denatured_rc ? denatured_rc[0] : null) : denatured_rc) ?: null
                tlist.collect { t -> [ ref, [ meta: gmeta, treated: t, untreated: untreated, denatured: denatured ] ] }
            }

        // Per reference: decide enablement and assemble positionally-paired control lists.
        def ch_nf_candidates = ch_nf_pairs
            .groupTuple()
            .map { ref, entries ->
                def base_meta     = entries[0].meta
                def treated_list  = entries.collect { entry -> entry.treated }
                def untreated_all = entries.collect { entry -> entry.untreated }
                def denatured_all = entries.collect { entry -> entry.denatured }
                def hasUntreated  = untreated_all.any { it }
                // rf-normfactor needs every treated paired with an untreated for Ding/Siegfried scoring.
                if (hasUntreated && untreated_all.any { !it }) {
                    def missing = entries.findAll { entry -> !entry.untreated }.collect { entry -> entry.meta.id }.sort().join(', ')
                    error("rf-normfactor for reference '${ref}': not every treated sample has a matched untreated control (missing for: ${missing}). Cross-experiment normalisation requires all-or-none untreated controls; add the missing controls or set --rfnorm_use_normfactor false.")
                }
                // AUTO: a reference is worth a shared factor once it has paired untreated controls or
                // more than one treated sample; an explicit 'true' forces it on regardless.
                def autoEnable     = hasUntreated || treated_list.size() > 1
                def enabled        = nfForceOn ? true : autoEnable
                // Denatured is optional, but must align 1:1 with treated to stay positionally paired; if
                // only some groups have a denatured control, drop it rather than mispair.
                def denatured_list = denatured_all.every { it } ? denatured_all : []
                def principle      = (base_meta.principle ?: '').toLowerCase()
                def scoringMethod  = principle == 'map'
                    ? (hasUntreated ? 3 : 4)
                    : (hasUntreated ? 1 : 2)
                def normMethod = resolveRfNormNormMethod(pipeline_config, scoringMethod)
                def nfmeta = base_meta + [
                    id                    : ref,
                    rfnorm_scoring_method : scoringMethod,
                    rfnorm_norm_method    : normMethod
                ]
                [ ref, enabled, nfmeta, treated_list, hasUntreated ? untreated_all : [], denatured_list ]
            }

        def ch_nf_input = ch_nf_candidates
            .filter { _ref, enabled, _meta, _t, _u, _d -> enabled }
            .map    { _ref, _enabled, nfmeta, t, u, d -> [ nfmeta, t, u, d, [] ] }

        RNAFRAMEWORK_RFNORMFACTOR(ch_nf_input)
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORMFACTOR.out.versions.first())

        // Complete per-reference factor map covering EVERY candidate reference, so all groups match
        // exactly one entry below. A reference gets a factor file only if it was enabled AND
        // rf-normfactor actually produced one; disabled references, and enabled ones where
        // rf-normfactor produced nothing (best-effort skip, e.g. coverage too low), get an empty slot
        // and fall back to per-sample normalisation.
        def ch_factor_by_ref = ch_nf_candidates
            .map { ref, _enabled, _meta, _t, _u, _d -> [ ref, true ] }
            .join(
                RNAFRAMEWORK_RFNORMFACTOR.out.factors.map { nfmeta, factors -> [ nfmeta.id.toString(), factors ] },
                remainder: true
            )
            .map { ref, _present, factors -> [ ref, factors ?: [] ] }

        ch_norm_input_final = ch_norm_input
            .map { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files ->
                [ resolveReferenceKey(gmeta, pipeline_config.organism), gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files ]
            }
            .combine(ch_factor_by_ref, by: 0)
            .map { _ref, gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, factors ->
                [ gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, factors ]
            }
    }

    def ch_xml
    def ch_plots
    def ch_rfnorm_log

    if (pipeline_config.rfnorm_chunk_size && !pipeline_config.transcriptome) {
        // Scatter: split the treated RC into transcript chunks; fan out a RFNORM job per chunk.
        // Only the first treated RC per group is split (groups with multiple treated samples are
        // uncommon; multi-treated support can be added if needed).
        // Pass both treated and untreated to SPLIT so matching chunks are extracted from both.
        // Transcript names in the chunk RC files must match between treated and untreated
        // for rf-norm to pair them correctly.
        //
        // Resolve the GTF for each group up front: SPLIT derives per-transcript lengths from it
        // (rf-rctools stats does not report lengths) to build the 4-column extraction BEDs.
        def ch_norm_input_with_gtf = ch_norm_input
            .combine(ch_reference_gtf_map)
            .map { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, gtf_map ->
                def reference_key = resolveReferenceKey(gmeta, pipeline_config.organism)
                def gtf_tuple     = gtf_map[reference_key]
                if (!gtf_tuple) {
                    error("No GTF resolved for reference '${reference_key}' (group '${gmeta.id}') — required to chunk RC files for rf-norm.")
                }
                [ gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, gtf_tuple[1] ]
            }

        def ch_norm_branched = ch_norm_input_with_gtf.multiMap { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, gtf ->
            for_split:    [ gmeta, treated_rcs instanceof List ? treated_rcs[0] : treated_rcs, untreated_rc ?: [], gtf ]
            for_controls: [ gmeta.id.toString(), denatured_rc, rci_files ]
        }

        RNAFRAMEWORK_RFRCTOOLS_SPLIT(ch_norm_branched.for_split)
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFRCTOOLS_SPLIT.out.versions.first())

        // Flatten treated chunk RC files; pair each with its matching untreated chunk by chunk index.
        def ch_treated_transposed = RNAFRAMEWORK_RFRCTOOLS_SPLIT.out.treated_chunks
            .transpose()
            .map { gmeta, chunk_rc ->
                def chunk_id    = (chunk_rc.name =~ /_chunk_(\d+)\.rc$/)[0][1]
                def chunk_gmeta = gmeta + [id: "${gmeta.id}_chunk_${chunk_id}".toString()]
                [ gmeta.id.toString(), chunk_id, chunk_gmeta, chunk_rc ]
            }

        def ch_untreated_transposed = RNAFRAMEWORK_RFRCTOOLS_SPLIT.out.untreated_chunks
            .transpose()
            .map { gmeta, chunk_rc ->
                def chunk_id = (chunk_rc.name =~ /_chunk_(\d+)_untreated\.rc$/)[0][1]
                [ gmeta.id.toString(), chunk_id, chunk_rc ]
            }

        def ch_chunked_inputs = ch_treated_transposed
            .join(ch_untreated_transposed, by: [0, 1], remainder: true)
            .combine(ch_norm_branched.for_controls, by: 0)
            .map { _group_key, _chunk_id, chunk_gmeta, treated_chunk, untreated_chunk, denatured_rc, rci_files ->
                // Chunking is mutually exclusive with --rfnorm_use_normfactor (guarded above), so the
                // factor slot is always empty on this path.
                [ chunk_gmeta, [treated_chunk], untreated_chunk ? [untreated_chunk] : [], denatured_rc ?: [], rci_files ?: [], [] ]
            }

        RNAFRAMEWORK_RFNORM(ch_chunked_inputs)
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())

        // Gather: re-collect outputs from all chunks back into one item per original group.
        ch_xml = RNAFRAMEWORK_RFNORM.out.xml
            .map { gmeta, xmls ->
                def orig_id = gmeta.id.replaceAll(/_chunk_\d+$/, '')
                [ orig_id, gmeta + [id: orig_id], xmls instanceof List ? xmls : [xmls] ]
            }
            .groupTuple(by: 0)
            .map { _orig_id, gmetas, xml_lists -> [ gmetas[0], xml_lists.flatten() ] }

        // Collect logs per original group — only the first chunk log is retained since
        // rf-norm's "covered" count is per-chunk; a proper aggregate would need summing.
        ch_rfnorm_log = RNAFRAMEWORK_RFNORM.out.log
            .map { gmeta, log ->
                def orig_id = gmeta.id.replaceAll(/_chunk_\d+$/, '')
                [ orig_id, gmeta + [id: orig_id], log ]
            }
            .groupTuple(by: 0)
            .map { _orig_id, gmetas, logs -> [ gmetas[0], logs[0] ] }

        ch_plots = RNAFRAMEWORK_RFNORM.out.plots
            .map { gmeta, plots ->
                def orig_id = gmeta.id.replaceAll(/_chunk_\d+$/, '')
                [ orig_id, gmeta + [id: orig_id], plots instanceof List ? plots : [plots] ]
            }
            .groupTuple(by: 0)
            .map { _orig_id, gmetas, plot_lists -> [ gmetas[0], plot_lists.flatten() ] }

    } else {
        RNAFRAMEWORK_RFNORM(ch_norm_input_final)
        ch_versions   = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
        ch_xml        = RNAFRAMEWORK_RFNORM.out.xml
        ch_plots      = RNAFRAMEWORK_RFNORM.out.plots
        ch_rfnorm_log = RNAFRAMEWORK_RFNORM.out.log
    }

    emit:
    xml         = ch_xml          // channel: [ val(meta), path(xml) ]
    plots       = ch_plots        // channel: [ val(meta), path(plots) ]
    rfnorm_log  = ch_rfnorm_log   // channel: [ val(meta), path(log) ]
    norm_groups = ch_rc_by_group  // channel: [ group, condition, val(meta), path(rc), path(rci) ]
    versions    = ch_versions     // channel: [ path(versions.yml) ]
}
