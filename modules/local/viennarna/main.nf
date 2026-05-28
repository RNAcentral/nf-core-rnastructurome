process VIENNARNA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.rnaframework_container

    input:
    tuple val(meta), path(fold_dir), path(xml, stageAs: "xml_input*/*"), path(drawn_ids)

    output:
    tuple val(meta), path("${prefix}_2D_structures/*.svg"), optional: true, emit: plots
    path "versions.yml",                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def rnaplot = task.ext.rnaplot ?: 'RNAplot'
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_2D_structures

    sort "${drawn_ids}" > .r2dt_drawn.txt

    for _db in ${fold_dir}/dotbracket/*.db; do
        [[ -f "\$_db" ]] || continue
        _id=\$(basename "\$_db" .db)

        if grep -qxF "\$_id" .r2dt_drawn.txt 2>/dev/null; then
            continue
        fi

        # Build per-transcript SHAPE file from XML reactivities
        python3 - "\${_id}" xml_input*/*.xml <<'PYEOF' > "\${_id}.shape" || true
import sys, xml.etree.ElementTree as ET, math

tid = sys.argv[1]
reps = []
for xml_path in sys.argv[2:]:
    try:
        root = ET.parse(xml_path).getroot()
    except Exception:
        continue
    for t in root.iter('transcript'):
        if (t.get('id') or t.findtext('id')) != tid:
            continue
        raw = t.findtext('reactivity') or t.findtext('values') or ''
        vals = []
        for v in raw.strip().split(','):
            v = v.strip()
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(float('nan'))
        if vals:
            reps.append(vals)
        break

if not reps:
    sys.exit(0)

length = max(len(r) for r in reps)
for i in range(length):
    col = [r[i] for r in reps if i < len(r) and not math.isnan(r[i])]
    v = sum(col) / len(col) if col else float('nan')
    out = -999 if math.isnan(v) or v < 0 else v
    print(f"{i+1}\t{out}")
PYEOF

        "${rnaplot}" --output-format=svg < "\$_db" || true

        if [[ -f "\${_id}_ss.svg" ]]; then
            mv "\${_id}_ss.svg" "${prefix}_2D_structures/\${_id}.svg"

            # Colour nucleotides by SHAPE reactivity if data available
            if [[ -s "\${_id}.shape" ]]; then
                python3 - "\${_id}.shape" "${prefix}_2D_structures/\${_id}.svg" <<'PYEOF2' || true
import sys, re, math, xml.etree.ElementTree as ET
shape_file, svg_file = sys.argv[1], sys.argv[2]
_SVG_NS = 'http://www.w3.org/2000/svg'

# Exact RNAframework color scheme (Interface/ViennaRNA.pm)
def bubble_colour(r):
    if r is None or math.isnan(r) or r < 0:
        return '#B1B3B6'   # grey      – no data
    if r <= 0.40:
        return '#000000'   # black     – low
    if r < 0.70:
        return '#FFCD2F'   # yellow    – medium
    return '#9A2322'       # dark red  – high

def text_colour(bg):
    return '#ffffff' if bg in ('#000000', '#9A2322') else '#000000'

reactivities = {}
with open(shape_file) as fh:
    for line in fh:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            try:
                reactivities[int(parts[0]) - 1] = float(parts[1])
            except ValueError:
                pass

ET.register_namespace('', _SVG_NS)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
tree = ET.parse(svg_file)
root = tree.getroot()

# Fix font and basepair stroke colour to match RNAframework output
for style_el in root.iter(f'{{{_SVG_NS}}}style'):
    if style_el.text:
        css = style_el.text
        css = css.replace('SansSerif', 'Arial, Helvetica, sans-serif')
        css = re.sub(r'(\.basepairs\b[^}]*?)stroke\s*:\s*red', r'\1stroke: #808080', css)
        style_el.text = css

parent_map = {child: parent for parent in root.iter() for child in parent}

texts = [e for e in root.iter(f'{{{_SVG_NS}}}text')
         if e.get('class', '').strip() == 'nucleotide']

for i, text_el in enumerate(texts):
    val = reactivities.get(i)
    fill = bubble_colour(val if val is not None else float('nan'))

    x = float(text_el.get('x', 0))
    y = float(text_el.get('y', 0))

    # Circle offset matches RNAframework: cx = x+4, cy = y-4, r=8, opacity=0.75
    circle = ET.Element(f'{{{_SVG_NS}}}circle')
    circle.set('cx', f'{x + 4:.3f}')
    circle.set('cy', f'{y - 4:.3f}')
    circle.set('r', '8')
    circle.set('fill', fill)
    circle.set('fill-opacity', '0.75')
    circle.set('stroke', fill)
    circle.set('stroke-width', '0.8')

    parent = parent_map.get(text_el)
    if parent is not None:
        parent.insert(list(parent).index(text_el), circle)

    # Text stays at original position; only fill colour changes
    text_el.set('style', f'fill: {text_colour(fill)};')

# Add reactivity legend to bottom-right (matching RNAframework layout)
svg_w = float(root.get('width', '452'))
svg_h = float(root.get('height', '452'))
sc = 0.75
tx = svg_w - 80 * sc - 5
ty = svg_h - 54 * sc - 10
legend = ET.SubElement(root, f'{{{_SVG_NS}}}g')
legend.set('transform', f'translate({tx:.1f},{ty:.1f}) scale({sc})')
for idx, (color, label) in enumerate([
    ('#9A2322', '0.7+'),
    ('#FFCD2F', '0.4-0.7'),
    ('#000000', '0-0.4'),
    ('#B1B3B6', 'No data'),
]):
    yo = idx * 13
    r = ET.SubElement(legend, f'{{{_SVG_NS}}}rect')
    r.set('x', '0'); r.set('y', str(yo)); r.set('width', '12'); r.set('height', '12')
    r.set('fill', color); r.set('stroke', '#000000'); r.set('stroke-width', '0.5')
    t = ET.SubElement(legend, f'{{{_SVG_NS}}}text')
    t.set('x', '15'); t.set('y', str(yo + 10))
    t.set('style', 'font-family: Arial, Helvetica, sans-serif; font-size: 11px; fill: #000000;')
    t.text = label

tree.write(svg_file, xml_declaration=True, encoding='unicode')
PYEOF2
            fi
        fi
    done

    printf '"%s":\\n    viennarna: %s\\n' \\
        "${task.process}" \\
        "\$("${rnaplot}" --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo 'unknown')" \\
        > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_2D_structures
    touch ${prefix}_2D_structures/stub_ENST00000000001.svg
    printf '"%s":\\n    viennarna: 2.6.4\\n' "${task.process}" > versions.yml
    """
}
