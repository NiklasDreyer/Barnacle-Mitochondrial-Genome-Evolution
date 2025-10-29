#!/bin/bash
#SBATCH --job-name=mafft_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=results/logs/mafft.%j.out
#SBATCH --error=results/logs/mafft.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/alignments results/logs

for gene_fa in results/per_gene_fastas/*.fa; do
  gene=$(basename "${gene_fa}" .fa)
  out="results/alignments/${gene}.afa"
  mafft --auto --thread 8 --maxiterate 1000 "${gene_fa}" > "${out}"
done

conda deactivate
echo "MAFFT alignments complete."
