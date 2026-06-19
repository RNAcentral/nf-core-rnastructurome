//
// NORMALISE_REACTIVITIES — group RC files by cell_line/replicate/condition, pair treated
// with untreated controls (exact match, or fuzzy base-token fallback when enabled), then
// run rf-norm to produce per-base reactivity XML files.
//

include { RNAFRAMEWORK_RFNORM } from '../../../modules/local/rnaframework/norm/main'
include {
    cellLineBaseToken
    resolveRfNormNormMethod
} from '../../../workflows/rnastructurome_functions.nf'

workflow NORMALISE_REACTIVITIES {

    take:
    ch_rfcount_rc   // channel: [ val(meta), path(rc) ]
    ch_rfcount_rci  // channel: [ val(meta), path(rci) ]
    pipeline_config // map

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

    // Key every sample by its normalisation group (cell_line_replicate) and condition.
    def ch_rc_by_group = ch_rc_with_rci
        .map { meta, rc, rci ->
            if (!meta.cell_line || !meta.replicate) {
                error("Missing cell_line or replicate for sample '${meta.id}'. rf-norm requires both columns to pair samples safely.")
            }
            def group     = "${meta.cell_line}_${meta.replicate}".toString()
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
    def ch_group_meta = ch_rc_by_group
        .map     { group, condition, meta, _rc, _rci -> [ group, [ condition: condition, meta: meta ] ] }
        .groupTuple()
        .map     { group, entries ->
            def conditions = entries.collect { entry -> entry.condition }.toSet()
            if (!conditions.contains('treated') && conditions.contains('untreated')) {
                def offending = entries
                    .findAll { entry -> entry.condition == 'untreated' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }.sort().join(', ')
                error("rf-norm requires a treated sample for cell_line/replicate group '${group}'. Invalid samples: ${offending}")
            }
            if (conditions.contains('denatured') && (!conditions.contains('treated') || !conditions.contains('untreated'))) {
                def offending = entries
                    .findAll { entry -> entry.condition == 'denatured' }
                    .collect { entry -> "${entry.meta.id} (${entry.condition})" }.sort().join(', ')
                error("rf-norm requires treated and untreated samples for denatured controls in cell_line/replicate group '${group}'. Invalid samples: ${offending}")
            }
            [ group, entries[0].meta ]
        }

    // Fuzzy untreated pairing fallback (enabled by default; disable with --fuzzy_untreated_pairing false).
    // When a treated group has no exact cell_line+replicate untreated match, the pipeline
    // falls back to an untreated sample sharing the same cell_line base token (the portion
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
                [ cellLineBaseToken(meta.cell_line.toString()), meta.replicate.toString(), group, rc ]
            }

        // Treated groups with no direct untreated match — candidates for fuzzy lookup.
        def ch_treated_no_untreated = ch_treated
            .join(ch_group_meta)
            .join(ch_untreated, remainder: true)
            .filter { _group, _treated_rcs, _base_meta, untreated_rc -> !untreated_rc }
            .map    { group, _treated_rcs, base_meta, _untreated_rc ->
                [ cellLineBaseToken(base_meta.cell_line.toString()), base_meta.replicate.toString(), group ]
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
                    error("Ambiguous untreated fallback for '${group}': multiple untreated groups share the same cell_line base and replicate: ${candidateGroups}. Use --fuzzy_untreated_pairing false to disable fuzzy matching.")
                }
                log.warn "No exact untreated match for '${group}' — falling back to '${candidates[0].untreated_group}' (shared cell_line base token at same replicate). Set --fuzzy_untreated_pairing false to require exact matches."
                [ group, candidates[0].rc ]
            }

        // Resolved untreated: direct exact match OR fuzzy fallback.
        ch_resolved_untreated = ch_untreated.mix(ch_fallback_untreated)
    } else {
        // Strict mode: only exact cell_line+replicate matches are used.
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

    RNAFRAMEWORK_RFNORM (
        ch_norm_input
    )
    ch_versions = ch_versions.mix(RNAFRAMEWORK_RFNORM.out.versions.first())

    emit:
    xml         = RNAFRAMEWORK_RFNORM.out.xml     // channel: [ val(meta), path(xml) ]
    plots       = RNAFRAMEWORK_RFNORM.out.plots   // channel: [ val(meta), path(plots) ]
    rfnorm_log  = RNAFRAMEWORK_RFNORM.out.log     // channel: [ val(meta), path(log) ]
    norm_groups = ch_rc_by_group                  // channel: [ group, condition, val(meta), path(rc), path(rci) ]
    versions    = ch_versions                     // channel: [ path(versions.yml) ]
}
