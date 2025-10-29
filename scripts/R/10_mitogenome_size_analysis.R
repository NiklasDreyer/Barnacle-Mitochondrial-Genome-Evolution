#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(readr)

opt_list <- list(
  make_option(c("-i","--input"), default="results/ref_db/mitogenome_sizes.csv", help="Input CSV"),
  make_option(c("-o","--outdir"), default="results/figures", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(opt$input, show_col_types = FALSE)

shapiro_test <- by(df$length_bp, df$infraclass, function(x) if(length(x) >=3) shapiro.test(x)$p.value else NA)
print("Shapiro test p-values per infraclass:")
print(shapiro_test)

test_result <- aov(length_bp ~ infraclass, data=df)
summary(test_result)

p <- ggplot(df, aes(x=infraclass, y=length_bp, fill=infraclass)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha=0.6) +
  theme_bw() + labs(x="Infraclass", y="Mitogenome length (bp)") +
  theme(legend.position = "none")

ggsave(filename=file.path(opt$outdir,"mitogenome_size_boxplot.png"), plot=p, width=6, height=4, dpi=300)

if(all(sapply(split(df$length_bp, df$infraclass), function(x) length(x) >= 3))){
  tuk <- TukeyHSD(test_result)
  write.csv(as.data.frame(tuk$infraclass), file.path(opt$outdir,"anova_tukey_results.csv"))
} else {
  library(rstatix)
  pairwise <- df %>%
    pairwise_wilcox_test(length_bp ~ infraclass, p.adjust.method = "BH")
  write.csv(pairwise, file.path(opt$outdir,"pairwise_wilcox_results.csv"), row.names = FALSE)
}

message("Mitogenome size analysis completed. Plots and tests in: ", opt$outdir)
