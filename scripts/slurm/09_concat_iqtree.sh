#!/bin/bash
#SBATCH --job-name=iqtree_concat
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=results/logs/iqtree.%j.out
#SBATCH --error=results/logs/iqtree.%j.err

set -euo pipefail
conda activate mito_env

mkdir -p results/iqtree results/logs

python - <<'PY'
from pathlib import Path
align_dir = Path("results/alignments/gblocks")
out_concat = Path("results/iqtree/concat.phy")
part_file = Path("results/iqtree/partitions.txt")
files = sorted(align_dir.glob("*.afa.gb"))
seqs = {}
order = []
for f in files:
    name = f.stem
    order.append(name)
    with open(f) as fh:
        header=None
        seq=""
        for line in fh:
            if line.startswith(">"):
                if header:
                    seqs.setdefault(header, []).append(seq)
                header=line[1:].strip()
                seq=""
            else:
                seq+=line.strip()
        if header:
            seqs.setdefault(header, []).append(seq)
samples = sorted(seqs.keys())
concat = {s: "".join(seqs[s]) for s in samples}
with open(out_concat, "w") as out:
    out.write(f"{len(samples)} {len(next(iter(concat.values())))}\n")
    for s in samples:
        out.write(f"{s} {concat[s]}\n")
pos=1
with open(part_file, "w") as p:
    for name in order:
        length=len(next(iter(seqs.values()))[0])
        end=pos+length-1
        p.write(f"DNA, {name} = {pos}-{end}\n")
        pos=end+1
print("Wrote concatenated alignment and partition file.")
PY

iqtree2 -s results/iqtree/concat.phy -spp results/iqtree/partitions.txt -m MFP+MERGE -bb 10000 -alrt 10000 -bnni -T AUTO -pre results/iqtree/mito_phylogeny

conda deactivate
echo "IQ-TREE run complete. Check results/iqtree directory."
