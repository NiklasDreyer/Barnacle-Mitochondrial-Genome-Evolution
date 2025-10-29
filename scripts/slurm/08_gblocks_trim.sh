#!/bin/bash
#SBATCH --job-name=gblocks_trim
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/gblocks.%j.out
#SBATCH --error=results/logs/gblocks.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/alignments/gblocks results/logs

for aln in results/alignments/*.afa; do
  name=$(basename "${aln}" .afa)
  Gblocks "${aln}" -t=c -b4=5 -b5=h
  if [[ -f "${aln}-gb" ]]; then
    mv "${aln}-gb" "results/alignments/gblocks/${name}.afa.gb"
  fi
done

conda deactivate
echo "Gblocks trimming finished."
