#!/usr/bin/env Rscript

# Full end-to-end DESeq2 + annotation + enrichment pipeline
# Usage:
#   Rscript run_full_end2end.R <salmon_quant_dir> <results_dir>
#
# Assumptions:
# - Salmon outputs are in <salmon_quant_dir>/<sample>/quant.sf
# - If sample_info.csv exists in project root, it must have columns: sample,group
# - If not, groups are inferred from sample names: UnT* -> Control, Q* -> TreatedQ, R* -> TreatedR

suppressPackageStartupMessages({
  # Install missing packages if needed
  required_cran <- c("dplyr","ggplot2","tibble","pheatmap","RColorBrewer")
  required_bioc  <- c("tximport","DESeq2","apeglm","edgeR","clusterProfiler","biomaRt","AnnotationDbi","org.Mm.eg.db")
  install_if_missing <- function(pkgs, bioc=FALSE){
    for(p in pkgs){
      if(!suppressWarnings(requireNamespace(p, quietly=TRUE))){
        if(!bioc){
          install.packages(p, repos="https://cloud.r-project.org")
        } else {
          if(!suppressWarnings(requireNamespace("BiocManager", quietly=TRUE))) install.packages("BiocManager")
          BiocManager::install(p, ask=FALSE, update=FALSE)
        }
      }
    }
  }
  install_if_missing(required_cran, bioc=FALSE)
  install_if_missing(required_bioc, bioc=TRUE)

  # Load libraries
  library(tximport); library(DESeq2); library(apeglm); library(edgeR)
  library(clusterProfiler); library(biomaRt); library(AnnotationDbi); library(org.Mm.eg.db)
  library(dplyr); library(ggplot2); library(tibble); library(pheatmap); library(RColorBrewer)
})

# -------- Arguments & dirs --------
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) stop("Usage: Rscript run_full_end2end.R <salmon_quant_dir> <results_dir>")
quant_dir <- args[1]
results_dir <- args[2]
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
plots_dir <- file.path(results_dir, "plots"); dir.create(plots_dir, showWarnings=FALSE)

# -------- Find Salmon quant files --------
quant_files <- list.files(quant_dir, pattern="quant.sf", recursive=TRUE, full.names=TRUE)
if(length(quant_files)==0) stop("No quant.sf files found in ", quant_dir)
sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

message("Found samples: ", paste(names(quant_files), collapse=", "))

# -------- Build sample table (use sample_info.csv if present) --------
# If sample_info.csv exists in parent of quant_dir (or working directory), use it; else infer.
project_root <- normalizePath(dirname(normalizePath(quant_dir)), mustWork=FALSE)
sample_info_paths <- c(file.path(project_root, "sample_info.csv"), file.path(getwd(), "sample_info.csv"))
sample_info_file <- sample_info_paths[file.exists(sample_info_paths)][1]
if(!is.na(sample_info_file)){
  message("Using sample_info.csv: ", sample_info_file)
  sample_table <- read.csv(sample_info_file, stringsAsFactors=FALSE)
  if(!all(c("sample","group") %in% colnames(sample_table))) stop("sample_info.csv must contain columns: sample,group")
  sample_table <- sample_table[match(names(quant_files), sample_table$sample), ]
} else {
  message("No sample_info.csv found — inferring groups from sample names.")
  infer_group <- function(s){
    if(grepl("^UnT", s, ignore.case=FALSE)) return("Control")
    if(grepl("^Q", s, ignore.case=FALSE)) return("TreatedQ")
    if(grepl("^R", s, ignore.case=FALSE)) return("TreatedR")
    # fallback: if contains "unt" or "ctrl"
    if(grepl("unt|ctrl", s, ignore.case=TRUE)) return("Control")
    if(grepl("^TQ|Q", s, ignore.case=TRUE)) return("TreatedQ")
    if(grepl("^TR|R", s, ignore.case=TRUE)) return("TreatedR")
    return("Unknown")
  }
  sample_table <- data.frame(sample = names(quant_files), group = sapply(names(quant_files), infer_group), stringsAsFactors = FALSE)
}
rownames(sample_table) <- sample_table$sample
message("Sample groups:\n"); print(sample_table)

# -------- tximport: optionally collapse transcripts->genes using biomaRt mapping --------
# Approach: inspect the first quant.sf to see if rows are transcript IDs (ENSMUST...) or gene IDs (ENSMUSG...), strip versions.
first_q <- quant_files[1]
first_df <- read.delim(first_q, header=TRUE, stringsAsFactors=FALSE)
ids_first <- first_df$Name
# strip version
ids_first_clean <- sub("\\..*$","", ids_first)
id_type <- if(all(grepl("^ENSMUST", ids_first_clean))) "TRANSCRIPT" else if(all(grepl("^ENSMUSG", ids_first_clean))) "GENE" else "OTHER"
message("Detected ID type in quant.sf: ", id_type)

tx2gene <- NULL
if(id_type == "TRANSCRIPT"){
  message("Building tx2gene mapping via biomaRt (this may take a moment)...")
  ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")
  transcripts <- unique(ids_first_clean)
  # get mapping for all transcripts present across samples (we'll query the union)
  all_ids_raw <- unique(unlist(lapply(quant_files, function(f) sub("\\..*$","", read.delim(f, header=TRUE, stringsAsFactors=FALSE)$Name))))
  mapping <- getBM(filters="ensembl_transcript_id",
                   attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                   values=all_ids_raw,
                   mart=ensembl)
  if(nrow(mapping)==0) stop("biomaRt mapping returned 0 rows. Check transcript IDs / internet / Ensembl version.")
  tx2gene <- unique(mapping[, c("ensembl_transcript_id","ensembl_gene_id")])
  colnames(tx2gene) <- c("tx","gene")
  message("tx2gene entries: ", nrow(tx2gene))
}

# Now run tximport (collapsing to genes if tx2gene exists)
message("Running tximport...")
if(!is.null(tx2gene)){
  # prepare named vector of quant file paths for tximport (Salmon expects quant.sf)
  txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE)
} else {
  txi <- tximport(quant_files, type="salmon", txOut=TRUE, ignoreTxVersion=TRUE)
}
message("tximport complete; matrix dims: ", paste(dim(txi$counts), collapse=" x "))

# -------- Create DESeq2 object and pre-filter low counts --------
sample_table <- sample_table[rownames(sample_table) %in% colnames(txi$counts), , drop=FALSE]
dds <- DESeqDataSetFromTximport(txi, colData=sample_table, design=~group)

# Filtering: use edgeR::filterByExpr if available, else rowSums >= 10
message("Filtering low-count genes...")
keep_edgeR <- tryCatch({
  keep <- filterByExpr(round(assay(txi$counts)), group=sample_table$group)
  TRUE
}, error=function(e) FALSE)
if(keep_edgeR){
  keep <- filterByExpr(round(assay(txi$counts)), group=sample_table$group)
  dds <- dds[keep,]
  message("Kept ", sum(keep), " genes after filterByExpr.")
} else {
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  message("Kept ", sum(keep), " genes after rowSums>=10 filter.")
}

# -------- Run DESeq and LFC shrinkage (apeglm) --------
dds$group <- factor(dds$group, levels=c("Control","TreatedQ","TreatedR"))
dds <- DESeq(dds)
message("DESeq finished.")

# Contrasts we'll test
contrasts <- list(
  Control_vs_TreatedQ = c("group","TreatedQ","Control"),
  Control_vs_TreatedR = c("group","TreatedR","Control"),
  TreatedQ_vs_TreatedR = c("group","TreatedQ","TreatedR")
)

# helper: annotation via biomaRt (strip versions first)
ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")

annotate_ids <- function(ids){
  ids_clean <- sub("\\..*$","", ids)
  # if gene-level (ENSMUSG) query gene; if transcript-level query transcript -> gene -> symbol & entrez
  if(all(grepl("^ENSMUSG", ids_clean))){
    map <- getBM(filters="ensembl_gene_id",
                 attributes=c("ensembl_gene_id","mgi_symbol","entrezgene_id"),
                 values=ids_clean, mart=ensembl)
    colnames(map) <- c("ensembl_gene_id","symbol","entrez")
    return(map)
  } else if(all(grepl("^ENSMUST", ids_clean))){
    map <- getBM(filters="ensembl_transcript_id",
                 attributes=c("ensembl_transcript_id","ensembl_gene_id","mgi_symbol","entrezgene_id"),
                 values=ids_clean, mart=ensembl)
    colnames(map) <- c("transcript_id","ensembl_gene_id","symbol","entrez")
    # return mapping with transcript_id as key
    return(map)
  } else {
    # assume SYMBOL
    map <- getBM(filters="mgi_symbol",
                 attributes=c("ensembl_gene_id","mgi_symbol","entrezgene_id"),
                 values=ids_clean, mart=ensembl)
    colnames(map) <- c("ensembl_gene_id","symbol","entrez")
    return(map)
  }
}

# output containers
all_sig_lists <- list()

# For each contrast: LFC shrink, save all + sig, volcano with legend, norm counts for sig, boxplots for top hits
for(cname in names(contrasts)){
  ct <- contrasts[[cname]]
  message("Processing contrast: ", cname)
  # lfcShrink returns shrunken LFCs; use type="apeglm"
  res_shr <- tryCatch({
    lfcShrink(dds, contrast=ct, type="apeglm")
  }, error=function(e){
    # fallback to results() if apeglm fails
    message("apeglm shrink failed: ", e$message, " — using results()")
    results(dds, contrast=ct)
  })
  res_df <- as.data.frame(res_shr)
  res_df <- tibble::rownames_to_column(res_df, var="feature")  # feature = gene or transcript
  # annotate
  ann <- annotate_ids(res_df$feature)
  # Merge: if mapping has transcript_id join on transcript, else try gene
  if("transcript_id" %in% colnames(ann)){
    res_df$transcript_id <- sub("\\..*$","", res_df$feature)
    res_df <- left_join(res_df, ann, by=c("transcript_id"="transcript_id"))
  } else if("ensembl_gene_id" %in% colnames(ann)){
    res_df$ensembl_gene_id <- sub("\\..*$","", res_df$feature)
    res_df <- left_join(res_df, ann, by=c("ensembl_gene_id"="ensembl_gene_id"))
  } else {
    # no annotation columns -> just keep feature as id
    res_df$symbol <- NA; res_df$entrez <- NA
  }
  # Save all DEGs
  all_file <- file.path(results_dir, paste0(cname,"_all_DEGs.csv"))
  write.csv(res_df, all_file, row.names=FALSE)
  # define significant DEGs by padj, abs(log2FC), baseMean
  PADJ_CUTOFF <- 0.05
  LFC_CUTOFF <- 1.0   # user-adjustable, 1 = two-fold
  BASEMEAN_CUTOFF <- 10
  sig_df <- res_df %>% filter(!is.na(padj) & padj < PADJ_CUTOFF & !is.na(log2FoldChange) &
                                abs(log2FoldChange) >= LFC_CUTOFF & baseMean >= BASEMEAN_CUTOFF)
  sig_file <- file.path(results_dir, paste0(cname,"_sig_DEGs.csv"))
  write.csv(sig_df, sig_file, row.names=FALSE)
  message("Saved ", nrow(sig_df), " significant DEGs to ", sig_file)
  all_sig_lists[[cname]] <- sig_df

  # Volcano plot: color by significance, legend counting only significant genes
  res_plot <- res_df %>% mutate(
    sig_status = case_when(
      is.na(padj) ~ "ns",
      padj < PADJ_CUTOFF & log2FoldChange > 0 ~ "up",
      padj < PADJ_CUTOFF & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    )
  )
  sig_only <- res_plot %>% filter(sig_status %in% c("up","down"))
  total_sig <- nrow(sig_only); upN <- sum(sig_only$sig_status=="up"); downN <- sum(sig_only$sig_status=="down")
  legend_text <- paste0("Significant genes: ", total_sig, "\nUp: ", upN, "\nDown: ", downN)
  p <- ggplot(res_plot, aes(x=log2FoldChange, y=-log10(padj), color=sig_status)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c("up"="red","down"="green","ns"="grey")) +
    theme_bw() + xlab("log2 Fold Change") + ylab("-log10(adj p-value)") +
    ggtitle(paste0("Volcano: ", cname)) +
    theme(legend.position="right") +
    guides(color=guide_legend(title=legend_text))
  ggsave(filename=file.path(plots_dir, paste0(cname,"_volcano.png")), plot=p, width=7, height=5)
  
  # normalized counts (vst) and boxplots for top hits
  vst_mat <- vst(dds, blind=FALSE)
  if(nrow(sig_df) > 0){
    top_hits <- head(sig_df$feature[order(-abs(sig_df$log2FoldChange))], 10)
    top_hits_clean <- sub("\\..*$","", top_hits)
    norm_counts <- assay(vst_mat)[top_hits_clean, , drop=FALSE]
    write.csv(as.data.frame(norm_counts), file=file.path(results_dir, paste0(cname,"_sig_norm_counts_top_hits.csv")), row.names=TRUE)
    # boxplots per gene
    for(g in rownames(norm_counts)){
      df_box <- data.frame(expr = norm_counts[g, ], sample = colnames(norm_counts))
      df_box$group <- sample_table[df_box$sample, "group"]
      pbox <- ggplot(df_box, aes(x=group, y=expr)) + geom_boxplot() + geom_jitter(width=0.2) +
        ggtitle(paste0(g, " (", cname, ")")) + ylab("VST expression")
      ggsave(filename=file.path(plots_dir, paste0(cname,"_",g,"_boxplot.png")), plot=pbox, width=5, height=4)
    }
  }
}

# -------- LRT for genes changing across any group --------
message("Running LRT (all-group) to find genes that differ across groups...")
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt) %>% tibble::rownames_to_column("feature") %>% arrange(padj)
write.csv(res_lrt_df, file=file.path(results_dir, "LRT_all_groups_all_DEGs.csv"), row.names=FALSE)
sig_lrt <- res_lrt_df %>% filter(!is.na(padj) & padj < 0.05 & !is.na(log2FoldChange))
write.csv(sig_lrt, file=file.path(results_dir, "LRT_all_groups_sig_DEGs.csv"), row.names=FALSE)
# LRT volcano
plot_volcano_lrt <- function(res_df, out_png, title="LRT Volcano"){
  res_plot <- res_df %>% mutate(
    sig_status = case_when(
      is.na(padj) ~ "ns",
      padj < 0.05 & log2FoldChange > 0 ~ "up",
      padj < 0.05 & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    )
  )
  sig_only <- res_plot %>% filter(sig_status %in% c("up","down"))
  legend_text <- paste0("Significant genes: ", nrow(sig_only), "\nUp: ", sum(sig_only$sig_status=="up"),
                        "\nDown: ", sum(sig_only$sig_status=="down"))
  p <- ggplot(res_plot, aes(x=log2FoldChange, y=-log10(padj), color=sig_status)) +
    geom_point(alpha=0.6) + scale_color_manual(values=c("up"="red","down"="green","ns"="grey")) +
    theme_bw() + xlab("log2 FC") + ylab("-log10 adj p") + ggtitle(title) +
    theme(legend.position="right") + guides(color=guide_legend(title=legend_text))
  ggsave(out_png, plot=p, width=7, height=5)
}
plot_volcano_lrt(res_lrt_df, file.path(plots_dir, "LRT_all_groups_volcano.png"), "LRT: all groups")

# -------- PCA (vst) --------
vst_mat <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vst_mat, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_pca <- ggplot(pcaData, aes(PC1, PC2, color=group)) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"%")) + ylab(paste0("PC2: ",percentVar[2],"%")) +
  theme_bw() + ggtitle("PCA (VST)")
ggsave(file.path(plots_dir,"PCA_all_samples.png"), plot=p_pca, width=6, height=5)

# -------- Top30 correlation heatmap using LRT significant genes --------
top30 <- sig_lrt %>% arrange(desc(abs(log2FoldChange))) %>% slice(1:30) %>% pull(feature)
top30_clean <- sub("\\..*$","", top30)
if(length(top30_clean) > 0){
  mat <- assay(vst_mat)[top30_clean, , drop=FALSE]
  mat_z <- t(scale(t(mat)))
  pheatmap(mat_z, annotation_col = sample_table["group"],
           color = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
           main="Top30 LRT DEGs (z-scored)", filename=file.path(plots_dir,"Top30_correlation_heatmap.png"),
           width=6, height=6)
}

# -------- Functional enrichment on significant lists only --------
message("Running functional enrichment (GO BP + KEGG) on significant DEG lists...")
# helper mapping for a sig_df: derive entrez IDs using biomaRt mapping (strip versions)
map_to_entrez <- function(ids){
  ids_clean <- sub("\\..*$","", ids)
  if(all(grepl("^ENSMUST", ids_clean))){
    m <- getBM(filters="ensembl_transcript_id", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl)
    entrez <- unique(m$entrezgene_id)
  } else if(all(grepl("^ENSMUSG", ids_clean))){
    m <- getBM(filters="ensembl_gene_id", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl)
    entrez <- unique(m$entrezgene_id)
  } else {
    m <- getBM(filters="mgi_symbol", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl)
    entrez <- unique(m$entrezgene_id)
  }
  entrez[!is.na(entrez)]
}

enrich_and_save <- function(sig_df, label){
  if(is.null(sig_df) || nrow(sig_df)==0){
    message("No sig genes for ", label, " — skipping enrichment.")
    return(NULL)
  }
  entrez <- map_to_entrez(sig_df$feature)
  entrez <- unique(entrez)
  if(length(entrez)==0){
    message("No mapped ENTREZ IDs for ", label, " — skipping enrichment.")
    return(NULL)
  }
  # GO BP
  ego <- enrichGO(gene=entrez, OrgDb=org.Mm.eg.db, keyType="ENTREZID", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
  # KEGG
  ekegg <- enrichKEGG(gene=entrez, organism='mmu', pAdjustMethod="BH", qvalueCutoff=0.05)
  # save
  write.csv(as.data.frame(ego), file=file.path(results_dir, paste0(label,"_GO_BP.csv")), row.names=FALSE)
  write.csv(as.data.frame(ekegg), file=file.path(results_dir, paste0(label,"_KEGG.csv")), row.names=FALSE)
  if(nrow(as.data.frame(ego))>0){
    png(file.path(plots_dir, paste0(label,"_GO_BP_dotplot.png")), width=900, height=700)
    print(dotplot(ego, showCategory=20) + ggtitle(paste0("GO BP Top20: ", label)))
    dev.off()
  } else message("No GO BP terms (q<0.05) for ", label)
  if(nrow(as.data.frame(ekegg))>0){
    png(file.path(plots_dir, paste0(label,"_KEGG_dotplot.png")), width=900, height=700)
    print(dotplot(ekegg, showCategory=20) + ggtitle(paste0("KEGG Top20: ", label)))
    dev.off()
  } else message("No KEGG terms (q<0.05) for ", label)
}

# run enrichment for each contrast sig list and LRT
for(name in names(all_sig_lists)){
  sigdf <- all_sig_lists[[name]]
  enrich_and_save(sigdf, name)
}
enrich_and_save(sig_lrt, "LRT_all_groups")

message("✅ Full pipeline finished. Results in: ", results_dir, " ; plots in: ", plots_dir)
