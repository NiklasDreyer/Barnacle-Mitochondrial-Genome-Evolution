#!/bin/bash
#SBATCH --job-name=velvet_assemble
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/velvet.%j.out
#SBATCH --error=results/logs/velvet.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/velvet results/logs

KMER=99
for fq in results/trimmed/*_R1.trimmed.fastq.gz; do
  sample=$(basename "${fq}" | sed 's/_R1.trimmed.fastq.gz//')
  r1="results/trimmed/${sample}_R1.trimmed.fastq.gz"
  r2="results/trimmed/${sample}_R2.trimmed.fastq.gz"
  outdir="results/velvet/${sample}"
  mkdir -p "${outdir}"
  velveth "${outdir}" ${KMER} -fastq.gz -shortPaired -separate "${r1}" "${r2}"
  velvetg "${outdir}" -exp_cov auto -cov_cutoff auto -ins_length 300
  cp "${outdir}/contigs.fa" "results/velvet/${sample}_contigs.fa"
done

conda deactivate
echo "Velvet assemblies completed."
