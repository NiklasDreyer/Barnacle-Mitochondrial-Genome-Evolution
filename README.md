# Barnacle_mitochondrial_genome_evolution

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/) [![R](https://img.shields.io/badge/R-4.3-orange)](https://www.r-project.org/) [![Status](https://img.shields.io/badge/status-draft-orange)]()

<img width="4292" height="2148" alt="image" src="https://github.com/user-attachments/assets/2052f0ea-efda-4da3-8a11-4dc6cec663af" />

Project overview
-----------------
Reproducible pipeline to assemble, annotate, and analyze mitochondrial genomes of Thecostraca (barnacles). Associated with the paper "Mitochondrial genome evolution and rearrangements in barnacles (Crustacea:Thecostraca) revealved by sampling of enigmatic lineages" by Dreyer et al. (Under review). Includes scripts to QC/trim reads, assemble contigs, identify mitochondrial contigs, annotate mitogenomes (MitoZ / MITOS), build reference databases, create per-gene alignments, and infer phylogenies with IQ-TREE2. Downstream comparative and per-gene rate/saturation analyses are implemented in R.

Key goals
- Generate circularized mitogenomes from Illumina paired-end reads.
- Provide consistent annotation (MitoZ / MITOS) and extract protein-coding genes.
- Build a reference mitogenome database (Cirripedia + Malacostraca outgroups).
- Produce codon-partitioned DNA and amino-acid matrices for phylogenetics.
- Estimate ML phylogenies with ModelFinder+MERGE and robust support (10k UFBoot + 10k SH-aLRT).
- Run per-gene evolutionary analyses (saturation, root-to-tip rates, clade comparisons) in R.

Repository layout
- `config.yaml` — project configuration (samples, paths, parameters)
- `environment.yml` — conda environment specification
- `LICENSE` — MIT license
- `.gitignore`
- `scripts/`  
  - `slurm/` — SLURM job scripts (QC/trim, assembly, BLAST, annotate, align, IQ-TREE)  
  - `bin/` — helper bash utilities (contig extraction, gene renaming)  
  - `R/` — R scripts for statistical analyses and plotting
- `data/` — place raw reads here (do not commit large raw data to git; use external storage)
- `results/` — pipeline outputs (assemblies, annotations, alignments, trees, figures)

Quick start (minimal)
1. Put raw read pairs in `data/raw_reads/` and update `config.yaml` so sample names and filenames match.
2. Create the environment:
   - Recommended: use `mamba` if available for faster dependency resolution:
     ```
     mamba env create -f environment.yml
     conda activate mito_env
     ```
3. Run QC + trim (cluster):
   ```
   sbatch scripts/slurm/00_qc_trim.slurm
   ```
4. Run subsequent pipeline steps in order (or use the Snakemake wrapper if you add one):
   - `01_assembly_velvet.slurm` → `02_blast_identify.slurm` → extract candidate mito contigs → `03_mitoz_annotate.slurm`
   - Build reference DB with `05_build_ref_db.slurm` (fill `results/ref_db/accessions.txt` first)
   - Prepare per-gene FASTAs, align with MAFFT, trim with Gblocks
   - Concatenate and run IQ-TREE: `09_concat_iqtree.slurm`

Short example (local test)
- Create a small test sample (or subsample reads) and run locally to validate environment and paths:
  ```
  # create test read dir and update config.yaml to point to test files
  mkdir -p data/raw_reads
  # run QC/trimming script interactively to verify tools
  bash scripts/slurm/00_qc_trim.slurm
  ```
- After annotation, extract protein-coding sequences for each taxon into `results/per_gene_fastas/` (one FASTA per gene, taxon names as headers), then:
  ```
  bash scripts/slurm/07_align_mafft.slurm
  bash scripts/slurm/08_gblocks_trim.slurm
  sbatch scripts/slurm/09_concat_iqtree.slurm
  ```

Expanded details / notes
- Trimming: Trimmomatic settings are in `config.yaml`. Make sure `TruSeq3-PE.fa` adapter file path is correct on your system or cluster.
- Assembly: Velvet k-mer and coverage settings can be tuned in `config.yaml` or the SLURM script. For mitochondrial genomes, a high k (e.g., 99) often works for 150 bp reads but test alternatives.
- Annotation: MitoZ is included in the conda environment (pip install). MITOS is typically used on the webserver — see `scripts/slurm/04_mitos_submit.sh` for guidance and where to store outputs.
- Per-gene FASTA extraction: The repo contains a template script (`scripts/bin/awk_rename_genes.sh`) but you will likely need to adapt it to MitoZ GFF/FASTA outputs. If you provide an example MitoZ output, the extraction script can be made robust.
- Phylogenetics: We use IQ-TREE2 with ModelFinder+MERGE; recommended settings are in `scripts/slurm/09_concat_iqtree.slurm`. Adjust CPU and memory to your cluster limits.
- R analyses: Scripts in `scripts/R/` implement mitogenome size comparisons and per-gene evolutionary tests (phangorn/ape-based). Use `Rscript` to run them non-interactively.

Outputs you should expect
- `results/velvet/*_contigs.fa` — assembly contigs
- `results/mitoz/<sample>/` — MitoZ annotation outputs (GFF, annotated FASTA)
- `results/per_gene_fastas/` — per-gene FASTAs (one file per gene)
- `results/alignments/` and `results/alignments/gblocks/` — alignments and trimmed alignments
- `results/iqtree/mito_phylogeny.*` — IQ-TREE outputs (.treefile, .log, .iqtree)
- `results/gene_analysis/` — saturation stats, per-taxon root-to-tip distance tables, diagnostic figures
- `results/figures/` — publication-ready PNGs

Contributing
- Use branches for features/bugfixes. Example:
  ```
  git checkout -b feature/mitoz-parser
  # make edits, add tests
  git commit -am "Add robust MitoZ GFF parser"
  git push -u origin feature/mitoz-parser
  ```
- Open a Pull Request describing the change and link related issues/pull requests.

Contact / citation
- Primary author: Niklas Dreyer, # Barnacle_mitochondrial_genome_evolution

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/) [![R](https://img.shields.io/badge/R-4.3-orange)](https://www.r-project.org/) [![Status](https://img.shields.io/badge/status-draft-orange)]()

Short description
-----------------
Reproducible pipeline to assemble, annotate, and analyze mitochondrial genomes of Thecostraca (barnacles). Includes scripts to QC/trim reads, assemble contigs, identify mitochondrial contigs, annotate mitogenomes (MitoZ / MITOS), build reference databases, create per-gene alignments, and infer phylogenies with IQ-TREE2. Downstream comparative and per-gene rate/saturation analyses are implemented in R.

Key goals
- Generate circularized mitogenomes from Illumina paired-end reads.
- Provide consistent annotation (MitoZ / MITOS) and extract protein-coding genes.
- Build a reference mitogenome database (Cirripedia + Malacostraca outgroups).
- Produce codon-partitioned DNA and amino-acid matrices for phylogenetics.
- Estimate ML phylogenies with ModelFinder+MERGE and robust support (10k UFBoot + 10k SH-aLRT).
- Run per-gene evolutionary analyses (saturation, root-to-tip rates, clade comparisons) in R.

Repository layout
- `config.yaml` — project configuration (samples, paths, parameters)
- `environment.yml` — conda environment specification
- `LICENSE` — MIT license
- `.gitignore`
- `scripts/`  
  - `slurm/` — SLURM job scripts (QC/trim, assembly, BLAST, annotate, align, IQ-TREE)  
  - `bin/` — helper bash utilities (contig extraction, gene renaming)  
  - `R/` — R scripts for statistical analyses and plotting
- `data/` — place raw reads here (do not commit large raw data to git; use external storage)
- `results/` — pipeline outputs (assemblies, annotations, alignments, trees, figures)

Quick start (minimal)
1. Put raw read pairs in `data/raw_reads/` and update `config.yaml` so sample names and filenames match.
2. Create the environment:
   - Recommended: use `mamba` if available for faster dependency resolution:
     ```
     mamba env create -f environment.yml
     conda activate mito_env
     ```
3. Run QC + trim (cluster):
   ```
   sbatch scripts/slurm/00_qc_trim.slurm
   ```
4. Run subsequent pipeline steps in order (or use the Snakemake wrapper if you add one):
   - `01_assembly_velvet.slurm` → `02_blast_identify.slurm` → extract candidate mito contigs → `03_mitoz_annotate.slurm`
   - Build reference DB with `05_build_ref_db.slurm` (fill `results/ref_db/accessions.txt` first)
   - Prepare per-gene FASTAs, align with MAFFT, trim with Gblocks
   - Concatenate and run IQ-TREE: `09_concat_iqtree.slurm`

Short example (local test)
- Create a small test sample (or subsample reads) and run locally to validate environment and paths:
  ```
  # create test read dir and update config.yaml to point to test files
  mkdir -p data/raw_reads
  # run QC/trimming script interactively to verify tools
  bash scripts/slurm/00_qc_trim.slurm
  ```
- After annotation, extract protein-coding sequences for each taxon into `results/per_gene_fastas/` (one FASTA per gene, taxon names as headers), then:
  ```
  bash scripts/slurm/07_align_mafft.slurm
  bash scripts/slurm/08_gblocks_trim.slurm
  sbatch scripts/slurm/09_concat_iqtree.slurm
  ```

Expanded details / notes
- Trimming: Trimmomatic settings are in `config.yaml`. Make sure `TruSeq3-PE.fa` adapter file path is correct on your system or cluster.
- Assembly: Velvet k-mer and coverage settings can be tuned in `config.yaml` or the SLURM script. For mitochondrial genomes, a high k (e.g., 99) often works for 150 bp reads but test alternatives.
- Annotation: MitoZ is included in the conda environment (pip install). MITOS is typically used on the webserver — see `scripts/slurm/04_mitos_submit.sh` for guidance and where to store outputs.
- Per-gene FASTA extraction: The repo contains a template script (`scripts/bin/awk_rename_genes.sh`) but you will likely need to adapt it to MitoZ GFF/FASTA outputs. If you provide an example MitoZ output, the extraction script can be made robust.
- Phylogenetics: We use IQ-TREE2 with ModelFinder+MERGE; recommended settings are in `scripts/slurm/09_concat_iqtree.slurm`. Adjust CPU and memory to your cluster limits.
- R analyses: Scripts in `scripts/R/` implement mitogenome size comparisons and per-gene evolutionary tests (phangorn/ape-based). Use `Rscript` to run them non-interactively.

Outputs you should expect
- `results/velvet/*_contigs.fa` — assembly contigs
- `results/mitoz/<sample>/` — MitoZ annotation outputs (GFF, annotated FASTA)
- `results/per_gene_fastas/` — per-gene FASTAs (one file per gene)
- `results/alignments/` and `results/alignments/gblocks/` — alignments and trimmed alignments
- `results/iqtree/mito_phylogeny.*` — IQ-TREE outputs (.treefile, .log, .iqtree)
- `results/gene_analysis/` — saturation stats, per-taxon root-to-tip distance tables, diagnostic figures
- `results/figures/` — publication-ready PNGs

Contributing
- Use branches for features/bugfixes. Example:
  ```
  git checkout -b feature/mitoz-parser
  # make edits, add tests
  git commit -am "Add robust MitoZ GFF parser"
  git push -u origin feature/mitoz-parser
  ```
- Open a Pull Request describing the change and link related issues/pull requests.

Contact / citation
- Primary author: Niklas Dreyer, # Barnacle_mitochondrial_genome_evolution

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/) [![R](https://img.shields.io/badge/R-4.3-orange)](https://www.r-project.org/) [![Status](https://img.shields.io/badge/status-draft-orange)]()

Short description
-----------------
Reproducible pipeline to assemble, annotate, and analyze mitochondrial genomes of Thecostraca (barnacles). Includes scripts to QC/trim reads, assemble contigs, identify mitochondrial contigs, annotate mitogenomes (MitoZ / MITOS), build reference databases, create per-gene alignments, and infer phylogenies with IQ-TREE2. Downstream comparative and per-gene rate/saturation analyses are implemented in R.

Key goals
- Generate circularized mitogenomes from Illumina paired-end reads.
- Provide consistent annotation (MitoZ / MITOS) and extract protein-coding genes.
- Build a reference mitogenome database (Cirripedia + Malacostraca outgroups).
- Produce codon-partitioned DNA and amino-acid matrices for phylogenetics.
- Estimate ML phylogenies with ModelFinder+MERGE and robust support (10k UFBoot + 10k SH-aLRT).
- Run per-gene evolutionary analyses (saturation, root-to-tip rates, clade comparisons) in R.

Repository layout
- `config.yaml` — project configuration (samples, paths, parameters)
- `environment.yml` — conda environment specification
- `LICENSE` — MIT license
- `.gitignore`
- `scripts/`  
  - `slurm/` — SLURM job scripts (QC/trim, assembly, BLAST, annotate, align, IQ-TREE)  
  - `bin/` — helper bash utilities (contig extraction, gene renaming)  
  - `R/` — R scripts for statistical analyses and plotting
- `data/` — place raw reads here (do not commit large raw data to git; use external storage)
- `results/` — pipeline outputs (assemblies, annotations, alignments, trees, figures)

Quick start (minimal)
1. Put raw read pairs in `data/raw_reads/` and update `config.yaml` so sample names and filenames match.
2. Create the environment:
   - Recommended: use `mamba` if available for faster dependency resolution:
     ```
     mamba env create -f environment.yml
     conda activate mito_env
     ```
3. Run QC + trim (cluster):
   ```
   sbatch scripts/slurm/00_qc_trim.slurm
   ```
4. Run subsequent pipeline steps in order (or use the Snakemake wrapper if you add one):
   - `01_assembly_velvet.slurm` → `02_blast_identify.slurm` → extract candidate mito contigs → `03_mitoz_annotate.slurm`
   - Build reference DB with `05_build_ref_db.slurm` (fill `results/ref_db/accessions.txt` first)
   - Prepare per-gene FASTAs, align with MAFFT, trim with Gblocks
   - Concatenate and run IQ-TREE: `09_concat_iqtree.slurm`

Short example (local test)
- Create a small test sample (or subsample reads) and run locally to validate environment and paths:
  ```
  # create test read dir and update config.yaml to point to test files
  mkdir -p data/raw_reads
  # run QC/trimming script interactively to verify tools
  bash scripts/slurm/00_qc_trim.slurm
  ```
- After annotation, extract protein-coding sequences for each taxon into `results/per_gene_fastas/` (one FASTA per gene, taxon names as headers), then:
  ```
  bash scripts/slurm/07_align_mafft.slurm
  bash scripts/slurm/08_gblocks_trim.slurm
  sbatch scripts/slurm/09_concat_iqtree.slurm
  ```

Expanded details / notes
- Trimming: Trimmomatic settings are in `config.yaml`. Make sure `TruSeq3-PE.fa` adapter file path is correct on your system or cluster.
- Assembly: Velvet k-mer and coverage settings can be tuned in `config.yaml` or the SLURM script. For mitochondrial genomes, a high k (e.g., 99) often works for 150 bp reads but test alternatives.
- Annotation: MitoZ is included in the conda environment (pip install). MITOS is typically used on the webserver — see `scripts/slurm/04_mitos_submit.sh` for guidance and where to store outputs.
- Per-gene FASTA extraction: The repo contains a template script (`scripts/bin/awk_rename_genes.sh`) but you will likely need to adapt it to MitoZ GFF/FASTA outputs. If you provide an example MitoZ output, the extraction script can be made robust.
- Phylogenetics: We use IQ-TREE2 with ModelFinder+MERGE; recommended settings are in `scripts/slurm/09_concat_iqtree.slurm`. Adjust CPU and memory to your cluster limits.
- R analyses: Scripts in `scripts/R/` implement mitogenome size comparisons and per-gene evolutionary tests (phangorn/ape-based). Use `Rscript` to run them non-interactively.

Outputs you should expect
- `results/velvet/*_contigs.fa` — assembly contigs
- `results/mitoz/<sample>/` — MitoZ annotation outputs (GFF, annotated FASTA)
- `results/per_gene_fastas/` — per-gene FASTAs (one file per gene)
- `results/alignments/` and `results/alignments/gblocks/` — alignments and trimmed alignments
- `results/iqtree/mito_phylogeny.*` — IQ-TREE outputs (.treefile, .log, .iqtree)
- `results/gene_analysis/` — saturation stats, per-taxon root-to-tip distance tables, diagnostic figures
- `results/figures/` — publication-ready PNGs

Contributing
- Use branches for features/bugfixes. Example:
  ```
  git checkout -b feature/mitoz-parser
  # make edits, add tests
  git commit -am "Add robust MitoZ GFF parser"
  git push -u origin feature/mitoz-parser
  ```
- Open a Pull Request describing the change and link related issues/pull requests.

Contact / citation
- Primary author: Niklas Dreyer, Méta-plateforme de génomique intégrée, Centre de recherche Azrieli du CHU Sainte-Justine, Montréal, QC, Canada, e-mail: niklasdreyerchen@gmail.com or niklas.dreyer.hsj@ssss.gouv.qc.ca
- If you use these scripts in a publication, please cite the corresponding papers and this repo.
Optional next steps I can prepare for you
- A short Snakemake file to orchestrate the pipeline (resumable, parallel).
- A robust MitoZ/GFF3 → per-gene FASTA extractor (I can write it if you upload an example MitoZ output).
- A small GitHub Action to run R script smoke tests on push (lightweight sanity checks).
- If you use these scripts in a publication, please cite the corresponding papers and this repo.

Acknowledgements
The sampling for this project was generously assisted by Pei-Chen "Vanessa" Tsai and Yao-Fong Tsao (Biodiversity Research Center, Taipei, Taiwan). We thank all personnel at the Marine Science Center, Green Island Marine Research Station (Academia Sinica) as well as the Sesoko Station (University of the Ryukyus). Niklas Dreyer greatly acknowledges support from postdoctoral advisor Prof. Gonzalo Giribet from the Museum of Comparative Zoology, Harvard University, USA, where the bioinformatic analyses for this project were carried out.
