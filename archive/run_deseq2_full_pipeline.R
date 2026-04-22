#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# -------- Arguments --------
args <- commandArgs(trailingOnly=TRUE)
quant_dir <- args[1]    # Salmon quant folder
results_dir <- args[2]  # Output folder
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)

# -------- Salmon quant files --------
files <- list.files(quant_dir, pattern="quant.sf", recursive=TRUE, full.names=TRUE)
names(files) <- basename(dirname(files))
if(length(files)==0) stop("No quant.sf files found in ", quant_dir)

txi <- tximport(files, type="salmon", txOut=TRUE)

# -------- Sample metadata --------
sample_table <- data.frame(
  sample = names(files),
  group = c("Control","Control","TreatedQ","TreatedQ","TreatedR","TreatedR")
)
rownames(sample_table) <- sample_table$sample
sample_table$group <- factor(sample_table$group, levels=c("Control","TreatedQ","TreatedR"))

# -------- DESeq2 analysis --------
dds <- DESeqDataSetFromTximport(txi, colData=sample_table, design=~group)
dds <- DESeq(dds)

# -------- Volcano plot function with significant-gene legend --------
plot_volcano <- function(res_df, plot_file, title="Volcano Plot") {
  # Identify significance for coloring
  res_df <- res_df %>%
    mutate(sig_status = case_when(
      is.na(padj) ~ "ns",
      padj < 0.05 & log2FoldChange > 0 ~ "up",
      padj < 0.05 & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    ))
  
  # Counts only for significant genes
  sig_genes_df <- res_df %>% filter(sig_status %in% c("up","down"))
  total_sig <- nrow(sig_genes_df)
  up_genes <- sum(sig_genes_df$sig_status=="up")
  down_genes <- sum(sig_genes_df$sig_status=="down")
  legend_text <- paste0("Significant genes: ", total_sig,
                        "\nUpregulated: ", up_genes,
                        "\nDownregulated: ", down_genes)
  
  p <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig_status)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("up"="red","down"="green","ns"="grey")) +
    theme_bw() +
    xlab("log2 Fold Change") +
    ylab("-log10 Adjusted P-value") +
    ggtitle(title) +
    theme(legend.position="right") +
    guides(color=guide_legend(title=legend_text))
  
  ggsave(plot_file, plot=p, width=6, height=5)
}

# -------- Pairwise contrasts --------
contrast_list <- list(
  Control_vs_TreatedQ = c("group","TreatedQ","Control"),
  Control_vs_TreatedR = c("group","TreatedR","Control"),
  TreatedQ_vs_TreatedR = c("group","TreatedQ","TreatedR")
)

for(cname in names(contrast_list)){
  contrast <- contrast_list[[cname]]
  message("Running contrast: ", cname)
  
  res <- results(dds, contrast=contrast)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("transcript_id") %>%
    arrange(padj)
  
  # Save all DEGs
  write.csv(res_df, file=file.path(results_dir, paste0(cname,"_all_DEGs.csv")), row.names=FALSE)
  
  # Save significant DEGs
  sig_res <- res_df %>% filter(!is.na(padj) & padj<0.05)
  write.csv(sig_res, file=file.path(results_dir, paste0(cname,"_sig_DEGs.csv")), row.names=FALSE)
  
  # Volcano plot with significant-gene legend
  plot_volcano(res_df, file.path(results_dir, paste0(cname,"_volcano.png")), title=paste0("Volcano: ", cname))
}

# -------- Likelihood ratio test (LRT) across all three groups --------
message("Running LRT for all-group differences...")
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt) %>%
  tibble::rownames_to_column("transcript_id") %>%
  arrange(padj)

# Save all DEGs and significant DEGs
write.csv(res_lrt_df, file=file.path(results_dir, "LRT_all_groups_all_DEGs.csv"), row.names=FALSE)
sig_lrt <- res_lrt_df %>% filter(!is.na(padj) & padj<0.05)
write.csv(sig_lrt, file=file.path(results_dir, "LRT_all_groups_sig_DEGs.csv"), row.names=FALSE)

# Volcano plot for LRT
plot_volcano(res_lrt_df, file.path(results_dir,"LRT_all_groups_volcano.png"), title="Volcano: All groups (LRT)")

# -------- PCA plot --------
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  ggtitle("PCA of all samples")

ggsave(file.path(results_dir, "PCA_all_samples.png"), plot=p_pca, width=6, height=5)

# -------- Top30 significant DEGs correlation heatmap --------
top30_genes <- res_lrt_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice(1:30) %>%
  pull(transcript_id)

rld_mat <- assay(rld)[top30_genes, ]
rld_mat_scaled <- t(scale(t(rld_mat)))  # Z-score per gene

pheatmap(rld_mat_scaled,
         annotation_col = sample_table["group"],
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         main = "Top30 Significant DEGs Correlation",
         filename = file.path(results_dir, "Top30_DEGs_correlation_heatmap.png"),
         width=6, height=6)

message("✅ Full DESeq2 pipeline complete. Results, volcano plots, PCA, and correlation heatmap are in ", results_dir)
