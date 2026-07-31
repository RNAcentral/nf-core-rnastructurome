// NORMALISE_REACTIVITIES — group RC files by sample_group/replicate/condition, pair treated with untreated
// controls (exact or fuzzy base-token fallback), then run rf-norm to produce per-base reactivity XML files.

include { RNAFRAMEWORK_RFNORM         } from '../../../modules/local/rnaframework/norm/main'
include { RNAFRAMEWORK_RFNORMFACTOR   } from '../../../modules/local/rnaframework/normfactor/main'
include {
    sampleGroupBaseToken
    resolveRfNormScoreMethod
    resolveRfNormNormMethod
    resolveReferenceKey
} from '../../../workflows/rnastructurome_functions.nf'

workflow NORMALISE_REACTIVITIES {

    take:
    ch_rfcount_rc        // channel: [ val(meta), path(rc) ]
    ch_rfcount_rci       // channel: [ val(meta), path(rci) ]
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

    // Enforce rf-norm pairing rules and select the representative meta per group, defined before the
    // fallback channels. The denatured-without-untreated check is deferred to ch_norm_input assembly.
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

    // Fuzzy untreated pairing fallback (default on; disable with --fuzzy_untreated_pairing false). When a
    // treated group has no exact untreated match, fall back to one sharing the same sample_group base token
    // (e.g. "MDA-MB-231" from "MDA-MB-231_MTX") at the same replicate — erroring if ambiguous.
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

        // Cross-product treated-without-untreated × available untreated, filtered on base+rep match,
        // then grouped by treated group to validate uniqueness before selecting the fallback.
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
            def scoringMethod = resolveRfNormScoreMethod(pipeline_config, principle, hasUntreated)
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

    // Reference-wide fallback: if a reference has exactly one untreated control, reuse it for every
    // treated group with none of its own (e.g. one shared control, several treated replicates).
    // Ambiguous cases (2+ distinct controls) fall through to the all-or-none check below.
    if (pipeline_config.fuzzy_untreated_pairing as Boolean) {
        ch_norm_input = ch_norm_input
            .map { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files ->
                def ref = resolveReferenceKey(gmeta, pipeline_config.organism)
                [ ref, [ meta: gmeta, treated: treated_rcs, untreated: untreated_rc, denatured: denatured_rc, rci: rci_files ] ]
            }
            .groupTuple()
            .flatMap { ref, entries ->
                // untreated may arrive as a bare Path or a single-element List (glob-typed process
                // output) — normalise to a bare Path (or null) before comparing/reusing.
                def untreatedFiles = entries.collect { entry ->
                    def u = entry.untreated
                    (u instanceof List ? (u ? u[0] : null) : u) ?: null
                }
                def distinct       = untreatedFiles.findAll { u -> u }.unique { u -> u.name }
                def missingIdx     = (0..<entries.size()).findAll { i -> !untreatedFiles[i] }
                if (distinct.size() == 1 && missingIdx) {
                    def missingIds = missingIdx.collect { i -> entries[i].meta.id }.sort()
                    log.warn "rf-norm reference '${ref}': only one untreated control ('${distinct[0].name}') is available — reusing it for treated group(s) with no untreated of their own: ${missingIds.join(', ')}. Set --fuzzy_untreated_pairing false to disable this fallback."
                }
                entries.withIndex().collect { entry, i ->
                    def gmeta        = entry.meta
                    def untreated_rc = untreatedFiles[i]
                    if (!untreated_rc && distinct.size() == 1) {
                        untreated_rc      = distinct[0]
                        def principle     = (gmeta.principle ?: '').toLowerCase()
                        def scoringMethod = resolveRfNormScoreMethod(pipeline_config, principle, true)
                        def normMethod    = resolveRfNormNormMethod(pipeline_config, scoringMethod)
                        gmeta = gmeta + [
                            rfnorm_has_untreated  : true,
                            rfnorm_scoring_method : scoringMethod,
                            rfnorm_norm_method    : normMethod
                        ]
                    }
                    [ gmeta, entry.treated, untreated_rc ?: [], entry.denatured, entry.rci ]
                }
            }
    }

    // Cross-experiment normalisation via rf-normfactor: derives one set of transcriptome-wide factors per
    // reference, fed to every group's rf-norm via -nf for a common scale (vs. per-sample box-plot).
    // Enablement is per reference: true=always on, false=always off, null=AUTO (on only when the reference
    // has >1 treated sample to cross-normalise). References that stay off normalise independently.
    def nfRaw      = pipeline_config.rfnorm_use_normfactor
    def nfForceOn  = (nfRaw != null) && (nfRaw.toString().toLowerCase() in ['true', '1', 'yes'])
    def nfForceOff = (nfRaw != null) && (nfRaw.toString().toLowerCase() in ['false', '0', 'no'])

    def ch_norm_input_final
    if (nfForceOff) {
        // rf-normfactor explicitly disabled.
        ch_norm_input_final = ch_norm_input
            .map { gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files ->
                [ gmeta, treated_rcs, untreated_rc, denatured_rc, rci_files, [] ]
            }
    } else {
        // Build per-reference candidates from ch_norm_input (already pairs treated with resolved
        // untreated/denatured). rf-normfactor pairs -t/-u/-d positionally, so emit one entry per
        // treated RC carrying its own controls, keeping the pairing intact regardless of groupTuple order.
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
                // groupTuple emits entries in run-varying arrival order, so sort by treated RC filename
                // to keep the staged lists stable — otherwise resume cache misses every time.
                def ordered       = entries.sort(false) { a, b -> a.treated.name <=> b.treated.name }
                def base_meta     = ordered[0].meta
                def treated_list  = ordered.collect { entry -> entry.treated }
                def untreated_all = ordered.collect { entry -> entry.untreated }
                def denatured_all = ordered.collect { entry -> entry.denatured }
                def hasUntreated  = untreated_all.any { u -> u }
                // rf-normfactor needs every treated paired with an untreated for Ding/Siegfried scoring.
                if (hasUntreated && untreated_all.any { u -> !u }) {
                    def missing = ordered.findAll { entry -> !entry.untreated }.collect { entry -> entry.meta.id }.sort().join(', ')
                    error("rf-normfactor for reference '${ref}': not every treated sample has a matched untreated control (missing for: ${missing}). Cross-experiment normalisation requires all-or-none untreated controls; add the missing controls or set --rfnorm_use_normfactor false.")
                }
                // AUTO: cross-experiment normalisation only makes sense with >1 treated sample to put on a
                // common scale. A single treated sample (even with an untreated control, which only affects
                // scoring) has nothing to cross-normalise, so per-sample normalisation suffices. Explicit
                // 'true' forces it on regardless.
                def autoEnable     = treated_list.size() > 1
                def enabled        = nfForceOn ? true : autoEnable
                // Denatured is optional, but must align 1:1 with treated to stay positionally paired; if
                // only some groups have a denatured control, drop it rather than mispair.
                def denatured_list = denatured_all.every { d -> d } ? denatured_all : []
                def principle      = (base_meta.principle ?: '').toLowerCase()
                def scoringMethod  = resolveRfNormScoreMethod(pipeline_config, principle, hasUntreated)
                def normMethod = resolveRfNormNormMethod(pipeline_config, scoringMethod)
                // Fuzzy pairing can resolve several treated samples to the SAME untreated/denatured control,
                // so the positional order (with repeats) needed for -t/-u/-d is carried as names in meta,
                // while the files themselves are staged deduplicated by name.
                def untreated_order = hasUntreated ? untreated_all : []
                def nfmeta = base_meta + [
                    id                    : ref,
                    rfnorm_scoring_method : scoringMethod,
                    rfnorm_norm_method    : normMethod,
                    nf_treated_order      : treated_list.collect    { f -> f.name },
                    nf_untreated_order    : untreated_order.collect  { f -> f.name },
                    nf_denatured_order    : denatured_list.collect   { f -> f.name }
                ]
                def treated_staged   = treated_list.unique    { f -> f.name }
                def untreated_staged = untreated_order.unique  { f -> f.name }
                def denatured_staged = denatured_list.unique   { f -> f.name }
                [ ref, enabled, nfmeta, treated_staged, untreated_staged, denatured_staged ]
            }

        def ch_nf_input = ch_nf_candidates
            .filter { _ref, enabled, _meta, _t, _u, _d -> enabled }
            .map    { _ref, _enabled, nfmeta, t, u, d -> [ nfmeta, t, u, d, [] ] }

        RNAFRAMEWORK_RFNORMFACTOR(ch_nf_input)
        ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORMFACTOR.out.versions.first())

        // Complete per-reference factor map covering EVERY candidate reference. A reference gets a factor
        // file only if enabled AND rf-normfactor produced one; otherwise it gets an empty slot and falls
        // back to per-sample normalisation (e.g. best-effort skip on low coverage).
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

    RNAFRAMEWORK_RFNORM(ch_norm_input_final)
    ch_versions   = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())
    def ch_xml        = RNAFRAMEWORK_RFNORM.out.xml
    def ch_plots      = RNAFRAMEWORK_RFNORM.out.plots
    def ch_rfnorm_log = RNAFRAMEWORK_RFNORM.out.log

    emit:
    xml         = ch_xml          // channel: [ val(meta), path(xml) ]
    plots       = ch_plots        // channel: [ val(meta), path(plots) ]
    rfnorm_log  = ch_rfnorm_log   // channel: [ val(meta), path(log) ]
    norm_groups = ch_rc_by_group  // channel: [ group, condition, val(meta), path(rc), path(rci) ]
    versions    = ch_versions     // channel: [ path(versions.yml) ]
}
