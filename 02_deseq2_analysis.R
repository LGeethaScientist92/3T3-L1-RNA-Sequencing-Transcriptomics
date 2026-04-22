#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(stringr)
})

project_dir <- "/mnt/c/RNASeqProject"
salmon_dir <- file.path(project_dir, "salmon_quant")
results_dir <- file.path(project_dir, "deseq2_results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------
# Build file list
# ------------------------------------------
files <- list.files(salmon_dir, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
sample_names <- basename(dirname(files))
names(files) <- sample_names

# ------------------------------------------
# Sample metadata
# ------------------------------------------
sample_table <- data.frame(
  sample = sample_names,
  group = ifelse(str_detect(sample_names, "^UnT"), "Control",
          ifelse(str_detect(sample_names, "^R"), "TreatedR",
          ifelse(str_detect(sample_names, "^Q"), "TreatedQ", NA)))
)

sample_table$group <- factor(sample_table$group, levels = c("Control", "TreatedR", "TreatedQ"))
rownames(sample_table) <- sample_table$sample

print(sample_table)

# ------------------------------------------
# tximport
# ------------------------------------------
txi <- tximport(files, type = "salmon", txOut = TRUE)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ group
)

dds <- DESeq(dds)
saveRDS(dds, file = file.path(results_dir, "dds.rds"))

# ------------------------------------------
# Export contrasts
# ------------------------------------------
res_control_vs_r <- results(dds, contrast = c("group", "TreatedR", "Control"))
res_control_vs_q <- results(dds, contrast = c("group", "TreatedQ", "Control"))
res_r_vs_q       <- results(dds, contrast = c("group", "TreatedR", "TreatedQ"))

write.csv(as.data.frame(res_control_vs_r), file.path(results_dir, "DEGs_TreatedR_vs_Control.csv"))
write.csv(as.data.frame(res_control_vs_q), file.path(results_dir, "DEGs_TreatedQ_vs_Control.csv"))
write.csv(as.data.frame(res_r_vs_q),       file.path(results_dir, "DEGs_TreatedR_vs_TreatedQ.csv"))

message("✅ DESeq2 analysis completed.")