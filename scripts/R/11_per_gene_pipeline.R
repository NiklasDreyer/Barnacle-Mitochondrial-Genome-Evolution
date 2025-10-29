#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(seqinr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(optparse)

opt_list <- list(
 make_option(c("-a","--alndir"), default="results/alignments", help="Directory with per-gene alignments (.afa or .fa)"),
 make_option(c("-o","--outdir"), default="results/gene_analysis", help="Output directory"),
 make_option(c("-m","--mapping"), default="results/ref_db/mitogenome_metadata.csv", help="Metadata mapping file (taxon -> clade)")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

aln_files <- list.files(opt$alndir, pattern="\\.afa$|\\.fa$|\\.fasta$", full.names = TRUE)

results_list <- list()
saturation_list <- list()
rates_table <- data.frame()

for(f in aln_files){
  gene <- tools::file_path_sans_ext(basename(f))
  cat("Processing gene:", gene, "\n")
  aln <- read.dna(f, format="fasta")
  if(nrow(aln) < 4) {
    cat("Skipping gene", gene, " (<4 taxa)\n")
    next
  }
  pd <- phyDat(as.character(aln), type="DNA")
  dm <- dist.ml(pd, model="GTR")
  tr <- NJ(dm)
  tr <- unroot(tr)
  fit <- pml(tr, data=pd)
  fit_opt <- optim.pml(fit, model="GTR", optGamma=TRUE, optInv=TRUE, optQ=TRUE, control = pml.control(trace=0))
  fit_tree <- fit_opt$tree
  if(!is.rooted(fit_tree)){
    fit_tree <- midpoint(fit_tree)
  }
  tip_depths <- node.depth.edgelength(fit_tree)[1:Ntip(fit_tree)]
  names(tip_depths) <- fit_tree$tip.label
  df_rates <- data.frame(taxon = names(tip_depths), gene = gene, root2tip = tip_depths, stringsAsFactors = FALSE)
  rates_table <- bind_rows(rates_table, df_rates)
  pdist <- as.matrix(dist.dna(aln, model="raw"))
  k80d <- as.matrix(dist.dna(aln, model="K80"))
  patristic <- cophenetic.phylo(fit_tree)
  get_upper <- function(m) {
    m[upper.tri(m)]
  }
  x1 <- get_upper(pdist)
  x2 <- get_upper(k80d)
  y  <- get_upper(patristic)
  lm_raw <- lm(y ~ x1)
  lm_k80 <- lm(y ~ x2)
  sat <- data.frame(gene=gene,
                    slope_raw = coef(lm_raw)[2], r2_raw = summary(lm_raw)$r.squared,
                    slope_k80 = coef(lm_k80)[2], r2_k80 = summary(lm_k80)$r.squared,
                    lnLik = fit_opt$logLik, AIC = AIC(fit_opt),
                    gamma = fit_opt$gamma, inv = fit_opt$inv)
  saturation_list[[gene]] <- sat
  p1 <- ggplot() +
    geom_point(aes(x=x1, y=y), alpha=0.5) +
    geom_smooth(aes(x=x1, y=y), method="lm", color="red") +
    labs(x="Uncorrected p-distance", y="ML patristic distance", title=paste("Saturation:", gene))
  ggsave(plot=p1, filename=file.path(opt$outdir, paste0(gene,"_sat_raw.png")), width=5, height=4)
  p2 <- ggplot() +
    geom_point(aes(x=x2, y=y), alpha=0.5) +
    geom_smooth(aes(x=x2, y=y), method="lm", color="red") +
    labs(x="K80 distance", y="ML patristic distance", title=paste("Saturation K80:", gene))
  ggsave(plot=p2, filename=file.path(opt$outdir, paste0(gene,"_sat_k80.png")), width=5, height=4)
}

write.csv(rates_table, file=file.path(opt$outdir, "per_taxon_root_to_tip.csv"), row.names=FALSE)
sat_df <- bind_rows(saturation_list)
write.csv(sat_df, file=file.path(opt$outdir, "saturation_statistics.csv"), row.names=FALSE)

if(file.exists(opt$mapping)){
  meta <- read.csv(opt$mapping, stringsAsFactors = FALSE)
  rates_mapped <- left_join(rates_table, meta, by=c("taxon"="taxon"))
  write.csv(rates_mapped, file=file.path(opt$outdir,"rates_mapped_clades.csv"), row.names = FALSE)
}

cat("Per-gene pipeline finished. Outputs in:", opt$outdir, "\n")
