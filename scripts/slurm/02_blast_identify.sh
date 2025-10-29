#!/bin/bash
#SBATCH --job-name=blast_identify
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/blast_identify.%j.out
#SBATCH --error=results/logs/blast_identify.%j.err

set -euo pipefail
conda activate mito_env

BLASTDB="/path/to/custom_mito_db.fasta"   # change to your DB path or "nt"
mkdir -p results/blast

for contig in results/velvet/*_contigs.fa; do
  sample=$(basename "${contig}" | sed 's/_contigs.fa//')
  out="results/blast/${sample}_blast.tsv"
  if [[ -f "${BLASTDB}" ]]; then
    blastn -query "${contig}" -db "${BLASTDB}" -outfmt '6 qseqid sseqid pident length evalue bitscore sseq' -evalue 1e-10 -num_threads 8 > "${out}"
  else
    blastn -query "${contig}" -db nt -remote -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue 1e-10 -num_threads 8 > "${out}"
  fi
done

conda deactivate
echo "BLAST identification finished."
