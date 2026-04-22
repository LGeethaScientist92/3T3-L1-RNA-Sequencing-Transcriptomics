#!/usr/bin/env Rscript
# run_full_pipeline_safe.R
# Full safe pipeline: tximport -> DESeq2 -> contrasts + LRT -> plots -> enrichment
# Hard-coded sample sheet (as provided).

suppressPackageStartupMessages({
  # ensure packages exist (install if needed) - optional, remove if you manage packages manually
  required_cran <- c("dplyr","ggplot2","tibble","pheatmap","RColorBrewer")
  required_bioc  <- c("tximport","DESeq2","apeglm","edgeR","clusterProfiler","biomaRt","AnnotationDbi","org.Mm.eg.db")
  for(p in required_cran) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  for(p in required_bioc) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

  # load libs
  library(tximport); library(DESeq2); library(apeglm); library(edgeR)
  library(clusterProfiler); library(biomaRt); library(AnnotationDbi); library(org.Mm.eg.db)
  library(dplyr); library(ggplot2); library(tibble); library(pheatmap); library(RColorBrewer)
})

# ---------- Arguments ----------
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) stop("Usage: Rscript run_full_pipeline_safe.R <salmon_quant_dir> <results_dir>")
quant_dir <- args[1]
results_dir <- args[2]
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
plots_dir <- file.path(results_dir, "plots"); dir.create(plots_dir, showWarnings=FALSE)

# ---------- Hard-coded sample sheet ----------
sample_table <- data.frame(
  sample = c("Q501_S13","Q502_S14","R501_S11","R502_S12","UnT1_S9","UnT2_S10"),
  group  = c("TreatmentQ","TreatmentQ","TreatmentR","TreatmentR","Control","Control"),
  stringsAsFactors = FALSE
)
rownames(sample_table) <- sample_table$sample
sample_table$group <- factor(sample_table$group, levels=c("Control","TreatmentQ","TreatmentR"))
message("Using sample table:")
print(sample_table)

# ---------- Find Salmon quant.sf files ----------
quant_files_all <- list.files(quant_dir, pattern="quant.sf", recursive=TRUE, full.names=TRUE)
if(length(quant_files_all)==0) stop("No quant.sf files found in ", quant_dir)
# map by parent dir name
sample_names_all <- basename(dirname(quant_files_all))
names(quant_files_all) <- sample_names_all
# keep only hard-coded samples (warn if missing)
common_samples <- intersect(names(quant_files_all), sample_table$sample)
if(length(common_samples) == 0) stop("No matching sample folders found in salmon_quant for the hard-coded sample names.")
if(length(common_samples) < nrow(sample_table)) warning("Some hard-coded samples missing; proceeding with available samples: ", paste(common_samples, collapse=", "))
quant_files <- quant_files_all[common_samples]
sample_table <- sample_table[common_samples, , drop=FALSE]
message("Using samples:", paste(names(quant_files), collapse=", "))

# ---------- detect ID type in quant.sf and prepare tx2gene if needed ----------
first_sf <- read.delim(quant_files[[1]], header=TRUE, stringsAsFactors=FALSE)
ids_first <- first_sf$Name
ids_first_clean <- sub("\\..*$","", ids_first)
id_type <- if(all(grepl("^ENSMUST", ids_first_clean))) "TRANSCRIPT" else if(all(grepl("^ENSMUSG", ids_first_clean))) "GENE" else "OTHER"
message("Detected ID type: ", id_type)

tx2gene <- NULL
ensembl <- NULL
if(id_type == "TRANSCRIPT"){
  message("Building tx2gene mapping with biomaRt (may take a minute)...")
  ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")
  # collect unique transcript IDs across available samples (strip versions)
  all_transcripts <- unique(unlist(lapply(quant_files, function(f) sub("\\..*$","", read.delim(f, header=TRUE, stringsAsFactors=FALSE)$Name))))
  mapping <- getBM(filters="ensembl_transcript_id",
                   attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                   values=all_transcripts,
                   mart=ensembl)
  if(nrow(mapping)==0) stop("biomaRt returned 0 rows for transcripts. Check IDs / internet / Ensembl version.")
  tx2gene <- unique(mapping[, c("ensembl_transcript_id","ensembl_gene_id")])
  colnames(tx2gene) <- c("tx","gene")
  message("tx2gene rows: ", nrow(tx2gene))
} else {
  # if not transcript-level we can still set ensembl for later mapping
  ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")
}

# ---------- tximport ----------
message("Running tximport...")
if(!is.null(tx2gene)){
  txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
} else {
  txi <- tximport(quant_files, type="salmon", txOut=TRUE, ignoreTxVersion=TRUE)
}
message("tximport done; dim counts: ", paste(dim(txi$counts), collapse=" x "))

# ---------- DESeq2 dataset ----------
# ensure column order matches sample_table
sample_table <- sample_table[colnames(txi$counts), , drop=FALSE]
dds <- DESeqDataSetFromTximport(txi, colData=sample_table, design=~group)

# ---------- filter low-count genes ----------
message("Filtering low-count features...")
keep <- NULL
try({
  keep <- filterByExpr(round(assay(txi$counts)), group=sample_table$group)
}, silent=TRUE)
if(is.null(keep) || length(keep)==0){
  keep <- rowSums(counts(dds)) >= 10
  message("Using fallback filter: rowSums >= 10")
} else {
  message("Using edgeR::filterByExpr")
}
dds <- dds[keep,]
message("Kept ", nrow(dds), " features after filtering")

# ---------- ensure design is valid ----------
dds$group <- droplevels(dds$group)
dds$group <- factor(dds$group, levels=c("Control","TreatmentQ","TreatmentR"))
if(nlevels(dds$group) < 2) stop("Less than 2 groups in design after filtering — cannot continue")
message("Group sample counts:"); print(table(colData(dds)$group))

# ---------- run DESeq ----------
message("Running DESeq() ...")
dds <- DESeq(dds)
message("DESeq finished.")

# ---------- set up contrasts and results container ----------
contrast_list <- list(
  Control_vs_TreatmentQ    = c("group","TreatmentQ","Control"),
  Control_vs_TreatmentR    = c("group","TreatmentR","Control"),
  TreatmentQ_vs_TreatmentR = c("group","TreatmentQ","TreatmentR")
)
results_list <- list()   # SAFEGUARD: always exists

# plotting helper: volcano with legend of significant genes only
plot_volcano <- function(res_df, out_png, title="Volcano", padj_cut=0.05){
  if(nrow(res_df) == 0) {
    message("No rows to plot for ", title); return(NULL)
  }
  # some padj may be NA -> protect when computing -log10
  res_plot <- res_df %>% mutate(padj_plot = ifelse(is.na(padj), 1, padj))
  res_plot <- res_plot %>% mutate(sig_status = case_when(
    is.na(padj) ~ "ns",
    padj < padj_cut & log2FoldChange > 0 ~ "up",
    padj < padj_cut & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))
  sig_only <- res_plot %>% filter(sig_status %in% c("up","down"))
  legend_text <- paste0("Significant genes: ", nrow(sig_only),
                        "\nUpregulated: ", sum(sig_only$sig_status=="up"),
                        "\nDownregulated: ", sum(sig_only$sig_status=="down"))
  p <- ggplot(res_plot, aes(x=log2FoldChange, y=-log10(padj_plot), color=sig_status)) +
    geom_point(alpha=0.6, size=1) +
    scale_color_manual(values=c("up"="red","down"="green","ns"="grey")) +
    theme_bw() + xlab("log2 Fold Change") + ylab("-log10(adj p)") +
    ggtitle(title) + theme(legend.position="right") +
    guides(color=guide_legend(title=legend_text))
  ggsave(out_png, plot=p, width=7, height=5)
}

# ---------- run contrasts safely ----------
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1.0
BASEMEAN_CUTOFF <- 10

for(name in names(contrast_list)){
  ct <- contrast_list[[name]]
  message("Running contrast: ", name)
  res_shr <- tryCatch({
    lfcShrink(dds, contrast=ct, type="apeglm")
  }, error=function(e){
    message("  lfcShrink failed for ", name, " (", e$message, ") - falling back to results()")
    tryCatch(results(dds, contrast=ct), error=function(e2){
      message("  results() also failed for ", name, ": ", e2$message); return(NULL)
    })
  })
  if(is.null(res_shr)){
    message("Skipping ", name, " because results unavailable.")
    next
  }
  res_df <- as.data.frame(res_shr) %>% tibble::rownames_to_column(var="feature")
  results_list[[name]] <- res_df
  # save all
  write.csv(res_df, file=file.path(results_dir, paste0(name,"_all_DEGs.csv")), row.names=FALSE)
  # filter significant
  sig_df <- res_df %>% filter(!is.na(padj) & padj < PADJ_CUTOFF & !is.na(log2FoldChange) &
                                abs(log2FoldChange) >= LFC_CUTOFF & baseMean >= BASEMEAN_CUTOFF)
  write.csv(sig_df, file=file.path(results_dir, paste0(name,"_sig_DEGs.csv")), row.names=FALSE)
  message("  Saved ", nrow(sig_df), " significant DEGs for ", name)
  # volcano (safe)
  plot_volcano(res_df, file.path(plots_dir, paste0(name,"_volcano.png")), title=paste0("Volcano: ", name), padj_cut=PADJ_CUTOFF)
}

# ---------- LRT (all-group) ----------
message("Running LRT (all-group) ...")
dds_lrt <- tryCatch({
  DESeq(dds, test="LRT", reduced=~1)
}, error=function(e){
  message("LRT DESeq failed: ", e$message); return(NULL)
})
if(!is.null(dds_lrt)){
  res_lrt <- results(dds_lrt)
  res_lrt_df <- as.data.frame(res_lrt) %>% tibble::rownames_to_column(var="feature") %>% arrange(padj)
  results_list[["LRT_all_groups"]] <- res_lrt_df
  write.csv(res_lrt_df, file=file.path(results_dir, "LRT_all_groups_all_DEGs.csv"), row.names=FALSE)
  sig_lrt <- res_lrt_df %>% filter(!is.na(padj) & padj < PADJ_CUTOFF)
  write.csv(sig_lrt, file=file.path(results_dir, "LRT_all_groups_sig_DEGs.csv"), row.names=FALSE)
  plot_volcano(res_lrt_df, file.path(plots_dir,"LRT_all_groups_volcano.png"), title="Volcano: LRT All Groups", padj_cut=PADJ_CUTOFF)
} else {
  message("Skipping LRT outputs due to earlier error.")
}

# ---------- Functional enrichment (GO BP + KEGG) on significant lists only ----------
message("Running functional enrichment on significant lists...")
# ensure ensembl mart exists
if(is.null(ensembl)) ensembl <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")

map_to_entrez <- function(ids){
  ids_clean <- sub("\\..*$","", ids)
  if(all(grepl("^ENSMUST", ids_clean))){
    m <- tryCatch(getBM(filters="ensembl_transcript_id", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl), error=function(e) data.frame(entrezgene_id=character(0)))
    entrez <- unique(m$entrezgene_id)
  } else if(all(grepl("^ENSMUSG", ids_clean))){
    m <- tryCatch(getBM(filters="ensembl_gene_id", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl), error=function(e) data.frame(entrezgene_id=character(0)))
    entrez <- unique(m$entrezgene_id)
  } else {
    m <- tryCatch(getBM(filters="mgi_symbol", attributes=c("entrezgene_id"), values=ids_clean, mart=ensembl), error=function(e) data.frame(entrezgene_id=character(0)))
    entrez <- unique(m$entrezgene_id)
  }
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

enrich_and_save <- function(sig_df, label){
  if(is.null(sig_df) || nrow(sig_df) == 0){
    message("No significant genes for ", label, "; skipping enrichment.")
    return(NULL)
  }
  entrez <- map_to_entrez(sig_df$feature)
  entrez <- unique(entrez)
  if(length(entrez) == 0){
    message("Mapping to ENTREZ returned 0 IDs for ", label, "; skipping enrichment.")
    return(NULL)
  }
  # GO BP
  ego <- tryCatch(enrichGO(gene=entrez, OrgDb=org.Mm.eg.db, keyType="ENTREZID", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE),
                  error=function(e) { message("enrichGO failed: ", e$message); return(NULL) })
  # KEGG
  ekegg <- tryCatch(enrichKEGG(gene=entrez, organism='mmu', pAdjustMethod="BH", qvalueCutoff=0.05), error=function(e){ message("enrichKEGG failed: ", e$message); return(NULL) })
  # Save tables
  if(!is.null(ego) && nrow(as.data.frame(ego))>0) write.csv(as.data.frame(ego), file=file.path(results_dir, paste0(label,"_GO_BP.csv")), row.names=FALSE)
  if(!is.null(ekegg) && nrow(as.data.frame(ekegg))>0) write.csv(as.data.frame(ekegg), file=file.path(results_dir, paste0(label,"_KEGG.csv")), row.names=FALSE)
  # Dotplots (top20) if results exist
  if(!is.null(ego) && nrow(as.data.frame(ego))>0){
    png(file.path(plots_dir, paste0(label,"_GO_BP_dotplot.png")), width=900, height=700)
    print(dotplot(ego, showCategory=20) + ggtitle(paste0("GO BP Top20: ", label)))
    dev.off()
  } else message("No GO BP terms (q<0.05) for ", label)
  if(!is.null(ekegg) && nrow(as.data.frame(ekegg))>0){
    png(file.path(plots_dir, paste0(label,"_KEGG_dotplot.png")), width=900, height=700)
    print(dotplot(ekegg, showCategory=20) + ggtitle(paste0("KEGG Top20: ", label)))
    dev.off()
  } else message("No KEGG terms (q<0.05) for ", label)
}

# run enrichment for each contrast's sig list and LRT
for(name in names(results_list)){
  res_df <- results_list[[name]]
  if(is.null(res_df) || nrow(res_df) == 0) next
  sig_df <- res_df %>% filter(!is.na(padj) & padj < PADJ_CUTOFF & !is.na(log2FoldChange) & abs(log2FoldChange) >= LFC_CUTOFF & baseMean >= BASEMEAN_CUTOFF)
  enrich_and_save(sig_df, name)
}
# LRT
if(exists("sig_lrt") && !is.null(sig_lrt)) enrich_and_save(sig_lrt, "LRT_all_groups")

message("All done. Results in: ", results_dir, " ; plots in: ", plots_dir)
