# MERGE_WIG

Local module that prepares BigWig inputs from per-transcript WIG files and RNAFramework XML outputs.

Behavior:

- derives `chrom.sizes` from transcript `id` and `length` attributes in the XML files
- merges per-transcript WIG files into a single `*.merged.wig`
- strips per-file `track` headers and emits one combined `track type=wiggle_0` header
