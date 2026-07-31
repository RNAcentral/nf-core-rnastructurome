/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnastructurome — pure helper functions
    Extracted from workflows/rnastructurome.nf (Phase 1). No channel/process logic;
    only value transforms, file parsers, arg renderers, and MultiQC table builders.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

def normaliseEnsemblSpecies(value) {
    value
        ?.toString()
        ?.trim()
        ?.toLowerCase()
        ?.replaceAll(/[^a-z0-9_]+/, '_')
        ?.replaceAll(/^_+|_+$/, '')
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

// Resolve how one sample's reference artifact (kind = 'fasta' | 'gtf') is obtained, returning
// [ reference_key, "<scheme>::<value>", original_organism ] (scheme: path:: | ensembl:: | ncbi:: | none::).
def resolveReferenceResolution(meta, cfg, kind) {
    def reference_key     = resolveReferenceKey(meta, cfg.organism)
    def original_organism = meta.organism?.toString() ?: reference_key
    def genome_entry      = cfg.genomes?.containsKey(reference_key) ? cfg.genomes[reference_key] : null

    // 1. Explicit local path (user-supplied or from params.genomes)
    if (kind == 'fasta') {
        def local_fasta = cfg.fasta
        if (local_fasta) {
            return [ reference_key, "path::${local_fasta.toString()}", original_organism ]
        }
        def transcript_fasta = genome_entry?.transcript_fasta ?: genome_entry?.transcriptome ?: genome_entry?.cdna
        if (transcript_fasta) {
            return [ reference_key, "path::${transcript_fasta.toString()}", original_organism ]
        }
    } else {
        if (cfg.gtf) {
            return [ reference_key, "path::${cfg.gtf.toString()}", original_organism ]
        }
        def gtf_path = genome_entry?.gtf
        if (gtf_path) {
            return [ reference_key, "path::${gtf_path.toString()}", original_organism ]
        }
    }

    // 2. Explicit Ensembl species (per-genome override or species map)
    def explicit_ensembl = genome_entry?.ensembl_species ?: cfg.ensembl_species_map?.get(reference_key)
    if (explicit_ensembl) {
        return [ reference_key, "ensembl::${explicit_ensembl.toLowerCase()}", original_organism ]
    }

    // 3. Pre-configured NCBI accessions. For FASTA this skips Ensembl (avoids spurious 404
    // log lines); for GTF the annotation is the synthetic NCBI_GTF, signalled by none::.
    def ncbi_accessions = genome_entry?.ncbi_accessions ?: cfg.ncbi_accessions_map?.get(reference_key)
    if (ncbi_accessions) {
        if (kind == 'fasta') {
            def acc_str = (ncbi_accessions instanceof List) ? ncbi_accessions.join(',') : ncbi_accessions.toString()
            return [ reference_key, "ncbi::${acc_str}", original_organism ]
        }
        return [ reference_key, "none::", original_organism ]
    }

    // 4. Default: Ensembl first. A 404 there triggers the NCBI fallback downstream.
    if (kind == 'fasta' && !reference_key) {
        error("No organism specified for sample '${meta.id}'. Provide --fasta, --organism, or set params.genomes.")
    }
    return [ reference_key, "ensembl::${reference_key}", original_organism ]
}

// Decide whether a run should default to the transcriptome (Bowtie) route: NCBI (bacteria/viral)
// references have no introns, so STAR offers nothing — auto-enable when every reference resolves to
// NCBI (samplesheets are single-organism-class by contract). Returns false on any parse error so the
// normal validation path in PIPELINE_INITIALISATION still runs and reports.
def allReferencesUseNcbiRoute(samplesheetPath, schemaPath, cfg) {
    if (!samplesheetPath) return false
    try {
        def rows = samplesheetToList(samplesheetPath.toString(), schemaPath.toString())
        if (!rows) return false
        def organisms = rows.collect { row ->
            def meta = (row instanceof List) ? row[0] : row
            (meta?.organism ?: cfg.organism)?.toString()?.trim()
        }
        if (organisms.any { org -> !org }) return false
        return organisms.every { org ->
            resolveReferenceResolution([ id: 'route-probe', organism: org ], cfg, 'fasta')[1].startsWith('ncbi::')
        }
    } catch (Exception _ignored) {
        return false
    }
}

// Collect a queue channel of [key, value] pairs into a single value channel holding
// a [key: value] map (last write wins on duplicate keys). Empty input yields [:].
def collectToMap(ch_keyed) {
    // .collect() on an empty channel emits nothing; ifEmpty keeps the result a real
    // (reusable value) map so downstream .combine() isn't silently emptied (e.g. no --gtf).
    ch_keyed
        .map { key, value -> [ (key): value ] }
        .collect()
        .map { entries -> entries.inject([:]) { acc, entry -> acc + entry } }
        .ifEmpty([:])
}

// Build STAR_ALIGN inputs for a set of trimmed reads: pair each sample with its reference's STAR index +
// GTF (combine by:0) and emit [ [meta,reads], [idx_meta,index], [gtf_meta,gtf], ignore_gtf ].
// ignore_gtf skips --sjdbGTFfile at align time only when no GTF was resolved for that reference.
def buildStarAlignInputs(ch_trimmed, ch_star_index, ch_gtf, cfg) {
    def ch_keyed = ch_trimmed
        .map { meta, reads -> [ resolveReferenceKey(meta, cfg.organism), meta, reads ] }
    def ch_idx_keyed = ch_star_index
        .map { meta, index -> [ meta.id.toString(), meta, index ] }
    def ch_gtf_keyed = ch_gtf
        .map { meta, gtf -> [ meta.id.toString(), meta, gtf ] }
    def ch_star_ref = ch_idx_keyed
        .join(ch_gtf_keyed, remainder: true)
        .map { combined ->
            def ref_key  = combined[0]
            def idx_meta = combined[1]
            def index    = combined[2]
            def gtf_meta = combined.size() > 3 ? combined[3] : null
            def gtf      = combined.size() > 4 ? combined[4] : null
            [ ref_key, idx_meta, index, gtf_meta ?: [id: 'no_gtf'], gtf ?: [], gtf_meta != null && gtf != null ]
        }
    ch_keyed
        .combine(ch_star_ref, by: 0)
        .map { _ref_key, sample_meta, reads, idx_meta, index, gtf_meta, gtf, has_gtf ->
            [ [sample_meta, reads], [idx_meta, index], [gtf_meta, gtf], !has_gtf ]
        }
}

// Collapse per-sample reference resolutions to one resolution per reference_key,
// erroring if a single reference resolved inconsistently across its samples.
def uniqueReferenceResolution(ch_requests, label) {
    ch_requests
        .groupTuple()
        .map { reference_key, resolutions, organisms ->
            def unique_resolutions = resolutions.unique()
            if (unique_resolutions.size() != 1) {
                error("Multiple ${label} resolutions were detected for reference '${reference_key}': ${unique_resolutions.join(', ')}")
            }
            [ reference_key, unique_resolutions[0], organisms[0] ]
        }
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
    // rf-count-genome's "second-strand" assigns transcript strand from read1's own mapped orientation
    // (read1=sense); "first-strand" uses read2's orientation (read2=sense, e.g. dUTP/TruSeq-directional).
    // RSeQC's "forward" fraction ("1++,1--,2+-,2-+") is read1=sense, so it maps to 'second' here; its
    // "reverse" fraction ("1+-,1-+,2++,2--", the common dUTP pattern) is read2=sense, mapping to 'first'.
    if (forward > 0.7) return 'second'
    if (reverse > 0.7) return 'first'
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

    // Hide aligner-specific params for the route not in use.
    // Transcriptome route uses Bowtie (RT-stop) / Bowtie2 (MaP); genome route uses STAR.
    if (params.transcriptome) {
        hiddenKeys.addAll(['star_map_sjdb_overhang', 'star_multimap_nmax'])
    } else {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('bowtie_') || key.startsWith('bowtie2_') })
    }

    // Hide genome-route-specific params for transcriptome runs, and vice versa
    if (params.transcriptome) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfcount_genome_') })
    } else {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfcount_map_') })
    }

    // Hide jackknife params when jackknife is not configured
    if (!params.jackknife_reference) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfjackknife_') || key == 'jackknife_reference' })
    }

    // Hide rfeval params when rfeval is not configured
    if (!params.rfeval_reference) {
        hiddenKeys.addAll(summaryParams.values()
            .findAll { section -> section instanceof Map }
            .collectMany { section -> section.keySet() as List }
            .findAll { key -> key.startsWith('rfeval_') || key == 'rfeval_reference' })
    }

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
    def args = []
    if (pipeline_config.bowtie_all as Boolean) {
        args << '-a'
    } else if (pipeline_config.bowtie_k != null) {
        args << "-k ${pipeline_config.bowtie_k as Integer}"
    }
    if ((pipeline_config.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${pipeline_config.bowtie_trim5 as Integer}"
    }
    if ((pipeline_config.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${pipeline_config.bowtie_trim3 as Integer}"
    }
    args << '-l 28'
    if (pipeline_config.bowtie_v != null) {
        args << "-v ${pipeline_config.bowtie_v as Integer}"
    } else {
        args << "-n ${pipeline_config.bowtie_n as Integer}"
    }
    if (!(pipeline_config.bowtie_all as Boolean) && pipeline_config.bowtie_k == null && pipeline_config.bowtie_max != null) {
        args << "-m ${pipeline_config.bowtie_max as Integer}"
        if ((pipeline_config.bowtie_max as Integer) > 1) {
            args << '-a'
        }
    }
    args << '--best'
    args << '--strata'
    args << "--chunkmbs ${pipeline_config.bowtie_chunkmbs as Integer}"
    args.join(' ').trim()
}

def renderBowtie2Args(pipeline_config) {
    def args = []
    def preset = pipeline_config.bowtie2_preset?.toString()?.trim() ?: ''
    if (pipeline_config.bowtie_all as Boolean) {
        args << '-a'
    } else if (!preset && pipeline_config.bowtie_k != null) {
        args << "-k ${pipeline_config.bowtie_k as Integer}"
    }
    if ((pipeline_config.bowtie_trim5 as Integer) > 0) {
        args << "--trim5 ${pipeline_config.bowtie_trim5 as Integer}"
    }
    if ((pipeline_config.bowtie_trim3 as Integer) > 0) {
        args << "--trim3 ${pipeline_config.bowtie_trim3 as Integer}"
    }
    if (preset) {
        args << preset
        if (pipeline_config.bowtie2_softclip as Boolean) {
            args << "--ma ${pipeline_config.bowtie2_ma as Integer}"
        }
    } else {
        args << '-L 22'
        if (pipeline_config.bowtie2_softclip as Boolean) {
            args << '--local'
            args << "--ma ${pipeline_config.bowtie2_ma as Integer}"
        }
    }
    args << "--mp ${pipeline_config.bowtie2_mp}"
    args << "--dpad ${pipeline_config.bowtie2_dpad as Integer}"
    args << "--rdg ${pipeline_config.bowtie2_rdg}"
    args << "--rfg ${pipeline_config.bowtie2_rfg}"
    if (pipeline_config.bowtie2_dovetail as Boolean) {
        args << '--dovetail'
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
    def scoringMethod = resolveRfNormScoreMethod(pipeline_config, principle, hasUntreated)
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
    if (pipeline_config.rfnorm_raw as Boolean) args << '--raw'
    if (pipeline_config.rfnorm_pseudocount != null) args << "--pseudocount ${pipeline_config.rfnorm_pseudocount}"
    if (pipeline_config.rfnorm_ignore_lower_than_untreated as Boolean) args << '--ignore-lower-than-untreated'
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

// Returns the base sample_group identifier (portion before the first underscore), e.g.
// "MDA-MB-231_MTX" → "MDA-MB-231". Used by fuzzy untreated-pairing to match a shared root.
def sampleGroupBaseToken(String sample_group) {
    sample_group.tokenize('_')[0]
}

def resolveRfNormScoreMethod(pipeline_config, principle, hasUntreated) {
    def defaultMethod = principle == 'map' ? (hasUntreated ? 3 : 4) : (hasUntreated ? 1 : 2)
    if (pipeline_config.rfnorm_score_method == null) {
        return defaultMethod
    }

    def requestedMethod = pipeline_config.rfnorm_score_method as Integer
    if (!(requestedMethod in [1, 2, 3, 4])) {
        error("Unsupported rf-norm scoring method '${pipeline_config.rfnorm_score_method}'. Expected one of: 1, 2, 3, 4.")
    }

    requestedMethod
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

// Normalise the assorted shapes Nextflow collect() can produce (single, nested, flattened, or Map) into
// a uniform list of [id, map] entries. With a non-null `label`, an unrecognised shape errors; else [].
def normaliseMqcRows(rows, label = null) {
    if (rows instanceof Map) {
        return rows.entrySet().collect { entry -> [entry.key, entry.value] }
    }
    if (rows instanceof List && rows.size() == 2 && rows[1] instanceof Map && !(rows[0] instanceof List)) {
        return [rows]
    }
    if (rows instanceof List && rows.every { row -> row instanceof List && row.size() == 2 && row[1] instanceof Map }) {
        return rows
    }
    if (rows instanceof List && rows.size() % 2 == 0 && rows.collate(2).every { pair -> pair.size() == 2 && pair[1] instanceof Map }) {
        return rows.collate(2)
    }
    if (label) {
        error("Unexpected ${label} row structure: ${rows?.getClass()?.name} -> ${rows}")
    }
    []
}

def countProgressionMultiqc(rows) {
    def orderedRows = normaliseMqcRows(rows, 'count progression').sort { a, b -> a[0] <=> b[0] }
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
    def orderedRows = normaliseMqcRows(rows).sort { a, b -> a[0] <=> b[0] }
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
            covered        : [title: 'RF-norm Good Coverage',  description: 'Transcripts passing RF-norm normalisation (sufficient coverage)', scale: 'Greens', format: '{:,.0f}']
        ]
    )
}

def rffoldStatsMultiqc(rows) {
    buildSimpleMultiqcTable(
        rows,
        'nf-core-rnastructurome-rffold-stats',
        'nf-core/rnastructurome RF-fold Statistics',
        'Folded transcripts per fold group, and transcripts dropped from the consensus fold because they were not covered in every replicate.',
        [
            folded   : [title: 'Folded Transcripts',   description: 'Transcripts successfully folded by rf-fold', scale: 'Purples', format: '{:,.0f}'],
            discarded: [title: 'Discarded Transcripts', description: 'Transcripts normalised in a replicate but excluded from the consensus fold (not present in all replicates)', scale: 'Reds', format: '{:,.0f}']
        ]
    )
}

// Parse an rf-correlate matrix.csv into the off-diagonal summary for the MultiQC reproducibility table
// (replicate count, mean/min pairwise correlation). Header: Sample,<label0>,...; rows: <label_i>,<corr_i0>,...
def parseRfcorrelateMatrix(matrixFile) {
    def lines = matrixFile.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return [ replicates: 0, mean_corr: 0, min_corr: 0 ]
    }
    def labels = lines[0].split(',').drop(1)
    def values = []
    lines.drop(1).eachWithIndex { line, i ->
        def parts = line.split(',')
        // Upper-triangle off-diagonal only: column j > row i. parts[0] is the row label.
        ((i + 1)..<labels.size()).each { j ->
            def raw = (j + 1) < parts.size() ? parts[j + 1].toString().trim() : ''
            if (raw && raw.toLowerCase() != 'nan') {
                values << (raw as Double)
            }
        }
    }
    if (!values) {
        return [ replicates: labels.size(), mean_corr: 0, min_corr: 0 ]
    }
    def mean = values.sum() / values.size()
    [
        replicates: labels.size(),
        mean_corr : (Math.round(mean * 1000) / 1000.0),
        min_corr  : (Math.round(values.min() * 1000) / 1000.0)
    ]
}

// rows: per sample_group [id, [replicates, mean_pearson, mean_spearman]], already merged across
// both correlation methods in the CORRELATE_REPLICATES subworkflow (flattened by .collect()).
def rfCorrelateMultiqc(rows) {
    def normalised = normaliseMqcRows(rows)
    if (!normalised) {
        return ''
    }
    buildSimpleMultiqcTable(
        normalised,
        'nf-core-rnastructurome-rfcorrelate',
        'nf-core/rnastructurome Replicate Correlation',
        'Pairwise reactivity-profile correlation between replicates from rf-correlate (per sample group). Higher is more reproducible.',
        [
            replicates    : [title: 'Replicates',      description: 'Number of replicates compared', scale: 'Blues',  format: '{:,.0f}'],
            mean_pearson  : [title: 'Mean Pearson',     description: 'Mean pairwise replicate Pearson correlation (reactivity-capped, overall, transcriptome-wide)', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}'],
            mean_spearman : [title: 'Mean Spearman',    description: 'Mean pairwise replicate Spearman correlation (overall, transcriptome-wide)', scale: 'RdYlGn', min: 0, max: 1, format: '{:,.3f}']
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
