#!/bin/bash
#SBATCH --job-name=build_ref_db
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/build_ref_db.%j.out
#SBATCH --error=results/logs/build_ref_db.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/ref_db

ACC_FILE="results/ref_db/accessions.txt"
OUTFA="results/ref_db/cirripedia_mitogenomes.fasta"

if [[ ! -s ${ACC_FILE} ]]; then
  cat > ${ACC_FILE} <<'EOF'
# Replace with real accession IDs (Table S1)
NC_012345.1
NC_023456.1
EOF
  echo "Please edit ${ACC_FILE} to include real accession IDs (one per line)."
fi

while read -r acc; do
  [[ -z "${acc}" || "${acc}" =~ ^# ]] && continue
  efetch -db nuccore -format fasta -id "${acc}" >> "${OUTFA}"
done < "${ACC_FILE}"

seqkit replace -p " .*" -r "" "${OUTFA}" -o "${OUTFA}.tmp"
mv "${OUTFA}.tmp" "${OUTFA}"

conda deactivate
echo "Reference database created at ${OUTFA}"
