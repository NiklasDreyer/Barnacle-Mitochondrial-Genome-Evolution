#!/bin/bash
#SBATCH --job-name=mitoz_annot
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=results/logs/mitoz.%j.out
#SBATCH --error=results/logs/mitoz.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/mitoz results/logs

for cand in results/mito_candidates/*.fa; do
  sample=$(basename "${cand}" | sed 's/.fa//')
  MitoZ all --genetic_code 5 --clade Arthropoda --thread_number 8 --outprefix "results/mitoz/${sample}" --fastafile "${cand}"
done

conda deactivate
echo "MitoZ annotations completed."
