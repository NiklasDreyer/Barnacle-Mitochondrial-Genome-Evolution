#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)
library(rstatix)
library(optparse)

opt_list <- list(
 make_option(c("-i","--rates"), default="results/gene_analysis/rates_mapped_clades.csv", help="Rates mapped csv"),
 make_option(c("-o","--outdir"), default="results/figures", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(opt$rates, show_col_types = FALSE)

counts <- df %>% group_by(gene,clade) %>% summarise(n=n())
valid_genes <- counts %>% group_by(gene) %>% summarise(min_n=min(n)) %>% filter(min_n >=3) %>% pull(gene)
df2 <- df %>% filter(gene %in% valid_genes)

res_list <- list()
for(g in unique(df2$gene)){
  tmp <- df2 %>% filter(gene==g)
  m <- lm(root2tip ~ clade, data=tmp)
  res <- resid(m)
  sw <- shapiro_test(as.data.frame(res))
  res_list[[g]] <- list(gene=g, normal_p=sw$p)
}
res_df <- bind_rows(lapply(res_list, as.data.frame))
write.csv(res_df, file=file.path(opt$outdir,"residuals_normality_per_gene.csv"), row.names=FALSE)

p <- ggplot(df2, aes(x=clade, y=root2tip, fill=clade)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, alpha=0.6) + facet_wrap(~gene, scales="free_y") +
  theme_bw() + theme(legend.position = "none")

ggsave(filename=file.path(opt$outdir,"root2tip_by_clade_per_gene.png"), plot=p, width=12, height=8, dpi=300)
cat("Rate visualizations saved to", opt$outdir, "\n")
