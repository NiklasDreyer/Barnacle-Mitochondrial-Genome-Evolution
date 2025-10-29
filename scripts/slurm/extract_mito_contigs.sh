#!/usr/bin/env bash
# Usage: extract_mito_contigs.sh blast_tsv contigs_fa outfa
set -euo pipefail
blast_tsv=$1
contigs_fa=$2
outfa=$3

cut -f1 "${blast_tsv}" | sort | uniq > /tmp/hit_ids.txt
seqkit grep -f /tmp/hit_ids.txt "${contigs_fa}" -o "${outfa}"

echo "Wrote candidate mitochondrial contigs to ${outfa}"
