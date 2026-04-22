#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# -------- Arguments --------
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  message("⚠️ No arguments provided. Using defaults...")
  dds_rds <- "/mnt/c/RNASeqProject/results/dds.rds"      # saved DESeq2 object
  results_dir <- "/mnt/c/RNASeqProject/results/QC"
} else {
  dds_rds <- args[1]
  results_dir <- args[2]
}
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)

# -------- Load DESeq2 object --------
if (!file.exists(dds_rds)) {
  stop("❌ File not found: ", dds_rds, 
       "\nPlease run DESeq2 first and save your DESeqDataSet as an .rds file.")
}
dds <- readRDS(dds_rds)

# -------- Variance stabilizing transform --------
message("Performing variance stabilizing transformation (VST)...")
vst_mat <- tryCatch({
  vst(dds, blind=TRUE)
}, error=function(e) {
  message("⚠️ VST failed: ", e$message)
  NULL
})

# -------- PCA plot --------
if (!is.null(vst_mat)) {
  message("Plotting PCA ...")
  pca_data <- plotPCA(vst_mat, intgroup="group", returnData=TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group, label=name)) +
    geom_point(size=4, alpha=0.8) +
    geom_text_repel(show.legend=FALSE, size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw(base_size=14) +
    ggtitle("PCA of Samples (VST)")
  
  ggsave(file.path(results_dir, "PCA_plot.png"), plot=p, width=6, height=5)
  ggsave(file.path(results_dir, "PCA_plot.pdf"), plot=p, width=6, height=5)
} else {
  message("❌ Skipping PCA plot (no VST).")
}

# -------- Sample correlation heatmap --------
if (!is.null(vst_mat)) {
  message("Plotting sample correlation heatmap ...")
  mat <- assay(vst_mat)
  cor_mat <- cor(mat, method="pearson")
  
  if (all(dim(cor_mat) > 1) && all(!is.na(cor_mat))) {
    pheatmap(cor_mat,
             annotation_col=as.data.frame(colData(dds)["group"]),
             main="Sample-to-Sample Correlation",
             color=colorRampPalette(brewer.pal(9,"Blues"))(255),
             border_color=NA,
             filename=file.path(results_dir, "sample_correlation_heatmap.png"),
             width=6, height=5)
  } else {
    message("❌ Correlation matrix is empty or invalid, skipping heatmap.")
  }
} else {
  message("❌ Skipping heatmap (no VST).")
}

message("✅ PCA + heatmap complete. Results saved in ", results_dir)
