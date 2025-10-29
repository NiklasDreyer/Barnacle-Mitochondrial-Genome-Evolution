#!/usr/bin/env bash
# Normalize gene names from MitoZ/MITOS outputs and assemble per-gene FASTAs.
input_fasta_dir="results/mitoz"
out_per_gene_dir="results/per_gene_fastas"
mkdir -p "${out_per_gene_dir}"

# This script is a template. You should parse GFF/GenBank outputs from MitoZ and extract CDS sequences.
# Example expects per-sample pcg FASTAs in results/mitoz/<sample>_pcgs/*.fa
for f in results/mitoz/*_pcgs/*.fa; do
  base=$(basename "${f}")
  gene=$(echo "${base}" | cut -d'_' -f1)
  sample=$(echo "${base}" | cut -d'_' -f2 | sed 's/.fa//')
  gene2=$(echo "${gene}" | awk '{print toupper($0)}')
  mkdir -p "${out_per_gene_dir}"
  awk -v sample="${sample}" -v gene="${gene2}" 'BEGIN{ORS=""} /^>/{print ">"sample"\n";next}{print $0"\n"}' "${f}" >> "${out_per_gene_dir}/${gene2}.fa"
done

echo "Per-gene FASTAs written to ${out_per_gene_dir}"
