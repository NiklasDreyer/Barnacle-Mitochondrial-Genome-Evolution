```markdown
# Barnacle_mitochondrial_genome_evolution

Repository overview
-------------------
This repository contains a reproducible pipeline used to assemble, annotate, and analyze mitochondrial genomes for Thecostraca (barnacles), including newly assembled Facetotecta (Type AE*) and Ascothoracida (Dendrogaster ludwigi) mitogenomes. The pipeline is organized as modular SLURM scripts (bash) for HPC runs and R scripts for downstream comparative analyses and visualization.

Main goals
- Assemble and circularize mitogenomes from Illumina paired-end reads.
- Annotate mitogenomes using MitoZ and MITOS (where applicable).
- Build a reference mitogenome database from publicly available Cirripedia + outgroups.
- Prepare codon-partitioned DNA and amino-acid matrices, run MAFFT and Gblocks.
- Estimate maximum likelihood phylogenies (IQ-TREE2) with model selection and rigorous support (UFBoot + SH-aLRT).
- Perform comparative mitogenome size analyses and per-gene substitution / saturation / rate comparisons in R.

Repository structure
- config.yaml               — project configuration (samples, paths, parameters)
- environment.yml           — conda environment to reproduce the pipeline
- LICENSE                   — license file (MIT)
- .gitignore                — patterns to ignore
- scripts/
  - slurm/                  — SLURM job scripts (QC, trim, assembly, blast, annotate, alignment, phylogeny)
  - bin/                    — helper bash scripts
  - R/                      — R scripts for statistics and plotting
- data/                     — (empty; place raw reads, assemblies, reference FASTA here)
- results/                  — (empty; pipeline outputs)

Quick start
1. Clone or create the repository on GitHub as `Barnacle_mitochondrial_genome_evolution`.
2. Create conda environment: `conda env create -f environment.yml` then `conda activate mito_env`.
3. Edit `config.yaml` with absolute/relative paths to your raw fastq files and sample names.
4. Submit SLURM jobs in order:
   - `scripts/slurm/00_qc_trim.slurm` (FastQC + Trimmomatic)
   - `scripts/slurm/01_assembly_velvet.slurm` (Velvet de novo assembly / contig generation)
   - `scripts/slurm/02_blast_identify.slurm` (identify mitochondrial contigs)
   - `scripts/slurm/03_mitoz_annotate.slurm` (annotation + coverage plots)
   - `scripts/slurm/05_build_ref_db.slurm` (download NCBI mitogenomes and prepare reference)
   - `scripts/slurm/07_align_mafft.slurm` (MAFFT per-gene alignments)
   - `scripts/slurm/08_gblocks_trim.slurm` (Gblocks trimming)
   - `scripts/slurm/09_concat_iqtree.slurm` (IQ-TREE phylogeny)
5. Run R analyses in `scripts/R/` for mitogenome sizes and per-gene analyses:
   - `Rscript scripts/R/10_mitogenome_size_analysis.R`
   - `Rscript scripts/R/11_per_gene_pipeline.R`

Notes and recommendations
- Adjust SLURM resources (CPUs, mem, time) to match your cluster.
- Update module/conda-loading lines in SLURM scripts to match your environment.
- MitoZ and MITOS: MitoZ can be run locally (included in environment.yml via pip). MITOS is often run via web server — `scripts/slurm/04_mitos_submit.sh` contains guidance.
- The pipeline expects short reads in paired-end FASTQ with names matching the `config.yaml` sample list.

Contact / authors
- Primary author: Niklas Dreyer (repository intended for use by the project team).
- Please cite the repository and underlying papers when using these scripts.
```
