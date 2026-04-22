#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(tibble)
})

# ===============================
# INPUT
# ===============================
dds_file <- "results/dds.rds"
outdir <- "results/plots"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

dds <- readRDS(dds_file)

# ===============================
# VST
# ===============================
vsd <- vst(dds, blind=TRUE)
mat <- assay(vsd)

# ===============================
# PCA
# ===============================
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=4) +
  geom_text_repel(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

ggsave(file.path(outdir, "PCA.png"), p, width=6, height=5)

# ===============================
# SAMPLE CORRELATION
# ===============================
cor_mat <- cor(mat)
pheatmap(cor_mat,
         annotation_col=as.data.frame(colData(dds)["group"]),
         color=colorRampPalette(brewer.pal(9,"Blues"))(100),
         filename=file.path(outdir, "Sample_correlation.png"))

# ===============================
# CONTRASTS
# ===============================
dds <- DESeq(dds)

# Control vs Treated
dds_ct <- dds
dds_ct$group2 <- ifelse(dds$group == "Control", "Control", "Treated")
dds_ct$group2 <- factor(dds_ct$group2)

dds_ct <- DESeq(dds_ct, design=~group2)
res_ct <- results(dds_ct, contrast=c("group2","Treated","Control"))

# R vs Q
dds_rq <- dds[, dds$group %in% c("TreatedR","TreatedQ")]
dds_rq$group <- droplevels(dds_rq$group)
dds_rq <- DESeq(dds_rq)
res_rq <- results(dds_rq, contrast=c("group","TreatedQ","TreatedR"))

# ===============================
# VOLCANO FUNCTION
# ===============================
plot_volcano <- function(res, title, filename){
  df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(sig = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ))

  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color=sig)) +
    geom_point(alpha=0.7) +
    scale_color_manual(values=c("red","blue","grey")) +
    theme_minimal() +
    ggtitle(title)

  ggsave(filename, p, width=6, height=5)
}

plot_volcano(res_ct, "Control vs Treated", file.path(outdir,"Volcano_CT.png"))
plot_volcano(res_rq, "R vs Q", file.path(outdir,"Volcano_RQ.png"))

# ===============================
# HEATMAP TOP 30 FUNCTION
# ===============================
plot_heatmap <- function(res, vsd, filename){
  res <- res[!is.na(res$padj),]
  res <- res[order(res$padj),]
  top30 <- head(res, 30)

  genes <- rownames(top30)
  mat <- assay(vsd)[genes,]
  mat <- t(scale(t(mat)))

  pheatmap(mat,
           annotation_col=as.data.frame(colData(dds)["group"]),
           color=colorRampPalette(c("navy","white","firebrick3"))(100),
           filename=filename)
}

plot_heatmap(res_ct, vsd, file.path(outdir,"Heatmap_CT.png"))
plot_heatmap(res_rq, vsd, file.path(outdir,"Heatmap_RQ.png"))

# ===============================
# KEGG ENRICHMENT
# ===============================
run_kegg <- function(res, prefix){
  res <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(padj < 0.05)

  gene_ids <- bitr(res$gene,
                   fromType="ENSEMBL",
                   toType="ENTREZID",
                   OrgDb=org.Mm.eg.db)

  ekegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                      organism="mmu")

  df <- as.data.frame(ekegg)
  write.csv(df, file.path(outdir, paste0(prefix,"_KEGG.csv")), row.names=FALSE)

  if(nrow(df)>0){
    df$negLog10 <- -log10(df$p.adjust)

    p <- ggplot(head(df,30),
                aes(negLog10, reorder(Description,negLog10))) +
      geom_point(aes(size=Count, color=negLog10)) +
      theme_minimal() +
      labs(title=paste(prefix,"KEGG"))

    ggsave(file.path(outdir,paste0(prefix,"_KEGG.png")), p, width=8, height=6)
  }
}

run_kegg(res_ct, "CT")
run_kegg(res_rq, "RQ")

message("✅ ALL plots generated")