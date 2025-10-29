#!/bin/bash
#SBATCH --job-name=qc_trim
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/qc_trim.%j.out
#SBATCH --error=results/logs/qc_trim.%j.err

set -euo pipefail

CONDDA_ENV="mito_env"

# load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${CONDDA_ENV}

mkdir -p results/fastqc results/trimmed results/logs

# Loop over samples in config.yaml using a small Python helper
python - <<'PY'
import yaml, os, subprocess
cfg = yaml.safe_load(open("config.yaml"))
rdir = cfg['raw_reads_dir']
for s in cfg['sample_table']:
    sample = s['sample']
    r1 = os.path.join(rdir, s['r1'])
    r2 = os.path.join(rdir, s['r2'])
    print("Processing", sample)
    # FastQC
    subprocess.run(["fastqc", "-o", "results/fastqc", "-t", "4", r1, r2], check=True)
    # Trimmomatic (adjust adapter path if necessary)
    trim_r1 = f"results/trimmed/{sample}_R1.trimmed.fastq.gz"
    trim_r2 = f"results/trimmed/{sample}_R2.trimmed.fastq.gz"
    trim_r1_un = f"results/trimmed/{sample}_R1.unpaired.fastq.gz"
    trim_r2_un = f"results/trimmed/{sample}_R2.unpaired.fastq.gz"
    cmd = [
        "trimmomatic", "PE", "-threads", "4", r1, r2,
        trim_r1, trim_r1_un, trim_r2, trim_r2_un,
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
        "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
    ]
    subprocess.run(cmd, check=True)
PY

conda deactivate
echo "QC and trimming finished."
