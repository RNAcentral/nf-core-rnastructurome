# FASTA_SORT

Local module that sorts FASTA records lexicographically by sequence identifier.

Inputs:

- reference FASTA, plain or gzip-compressed

Behavior:

- reads all FASTA records
- sorts by the first token in each header line
- writes a normalized plain-text `.sorted.fa`

This is used to keep transcript references deterministic before index building and downstream counting.
