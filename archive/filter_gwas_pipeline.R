# ====================================================
# GWAS Catalog pipeline for obesity DEGs + Enrichment + Volcano plots
# ====================================================

# --------------------------
# Install/load dependencies
# --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c("gwasrapidd", "dplyr", "readr", "biomaRt", 
              "clusterProfiler", "org.Mm.eg.db", "ReactomePA", "ggplot2", "enrichplot")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("biomaRt", "clusterProfiler", "org.Mm.eg.db", "ReactomePA", "enrichplot")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

library(gwasrapidd)
library(dplyr)
library(readr)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggplot2)
library(enrichplot)

# --------------------------
# 1. Load your DEG list
# --------------------------
deg_list <- read_csv("LRT_all_groups_sig_DEGs.csv")

# detect Ensembl gene ID column
if ("gene_id" %in% colnames(deg_list)) {
  genes_ensembl <- unique(na.omit(deg_list$gene_id))
} else {
  ens_col <- grep("ensembl|gene", colnames(deg_list), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(ens_col)) stop("❌ No Ensembl column found in input file")
  genes_ensembl <- unique(na.omit(deg_list[[ens_col]]))
  message("Using column: ", ens_col, " as Ensembl gene IDs")
}

cat("Loaded", length(genes_ensembl), "Ensembl IDs\n")

# --------------------------
# 2. Map Ensembl -> Mouse Symbols + Entrez IDs
# --------------------------
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = genes_ensembl,
  mart = ensembl_mouse
)

deg_mapped <- deg_list %>%
  left_join(mapping, by = c("gene_id" = "ensembl_gene_id"))

genes_symbols <- unique(na.omit(deg_mapped$mgi_symbol))
cat("Mapped to", length(genes_symbols), "mouse gene symbols\n")

# --------------------------
# 3. Query GWAS Catalog for obesity-related traits
# --------------------------
traits <- c("obesity", "body mass index", "waist hip ratio",
            "type 2 diabetes", "fasting glucose", "insulin")

studies <- get_studies(efo_trait = traits)
associations <- get_associations(study_id = studies@studies$study_id)

gwas_genes <- unique(associations@genes$gene_name)
cat("Retrieved", length(gwas_genes), "unique human GWAS obesity-related genes\n")

# --------------------------
# 3b. Map human GWAS genes -> mouse orthologs
# --------------------------

# Human mart
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map human gene symbols -> human Ensembl IDs
human_ids <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = gwas_genes,
  mart = ensembl_human
)

# Now map human Ensembl IDs -> mouse orthologs
human2mouse <- getLDS(
  attributes = c("hgnc_symbol"),      # human symbols
  filters = "hgnc_symbol",
  values = gwas_genes,
  mart = ensembl_human,
  attributesL = c("mgi_symbol"),      # mouse symbols
  martL = ensembl_mouse
)

gwas_mouse_genes <- unique(na.omit(human2mouse$MGI.symbol))
cat("Mapped GWAS obesity-related genes to", length(gwas_mouse_genes), "mouse ortholog symbols\n")


# --------------------------
# 4. Intersect DEGs with GWAS-mapped mouse genes
# --------------------------
obesity_related <- deg_mapped %>%
  filter(mgi_symbol %in% gwas_mouse_genes)

write_csv(obesity_related, "GWAS_Obesity_related_mouse_DEGs.csv")
cat("Found", nrow(obesity_related), "obesity-related DEGs. Saved to GWAS_Obesity_related_mouse_DEGs.csv\n")

# --------------------------
# 5. Pathway enrichment on obesity-related subset + plots
# --------------------------
obesity_entrez <- na.omit(unique(obesity_related$entrezgene_id))

if (length(obesity_entrez) > 0) {
  # Named vector for logFC (for coloring plots)
  geneList <- obesity_related$logFC
  names(geneList) <- obesity_related$entrezgene_id
  
  # KEGG enrichment
  kegg_res <- enrichKEGG(
    gene = obesity_entrez,
    organism = "mmu",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  write.csv(as.data.frame(kegg_res), "GWAS_KEGG_enrichment_obesity_DEGs.csv", row.names = FALSE)
  
  # Reactome enrichment
  reactome_res <- enrichPathway(
    gene = obesity_entrez,
    organism = "mouse",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  write.csv(as.data.frame(reactome_res), "GWAS_Reactome_enrichment_obesity_DEGs.csv", row.names = FALSE)
  
  # --------------------------
  # Enrichment plots with logFC coloring
  # --------------------------
  if (nrow(as.data.frame(kegg_res)) > 0) {
    png("GWAS_KEGG_cnetplot.png", width = 2500, height = 2000, res = 300)
    print(cnetplot(kegg_res, foldChange = geneList, showCategory = 10, colorEdge = TRUE))
    dev.off()
  }
  
  if (nrow(as.data.frame(reactome_res)) > 0) {
    png("GWAS_Reactome_cnetplot.png", width = 2500, height = 2000, res = 300)
    print(cnetplot(reactome_res, foldChange = geneList, showCategory = 10, colorEdge = TRUE))
    dev.off()
  }
  
  cat("Enrichment results saved: KEGG + Reactome CSVs and plots (PNG)\n")
  
} else {
  cat("No GWAS-matching DEGs found for enrichment.\n")
}

# --------------------------
# 6. Volcano plots
# --------------------------
if ("logFC" %in% colnames(obesity_related) & "pvalue" %in% colnames(obesity_related)) {
  obesity_related <- obesity_related %>%
    mutate(sig = case_when(
      logFC > 1 & pvalue < 0.05 ~ "Up",
      logFC < -1 & pvalue < 0.05 ~ "Down",
      TRUE ~ "NS"
    ))
  
  # ---- GWAS-only volcano with top 30 labels ----
  top30 <- obesity_related %>%
    arrange(pvalue) %>%
    slice_head(n = 30)
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
  
  png("Volcano_GWAS_Obesity_DEGs.png", width = 2200, height = 2000, res = 300)
  ggplot(obesity_related, aes(x = logFC, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_minimal(base_size = 14) +
    labs(title = "Volcano Plot (GWAS-linked Obesity DEGs)",
         x = "log2 Fold Change",
         y = "-log10(p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top30, aes(label = mgi_symbol),
                    size = 3.5, max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.3,
                    segment.color = "black") +
    theme(legend.position = "top")
  dev.off()
  
  cat("Volcano plot (GWAS subset) saved with top 30 DEGs labeled: Volcano_GWAS_Obesity_DEGs.png\n")
}

# ---- Full DEG volcano with GWAS DEGs highlighted ----
if ("logFC" %in% colnames(deg_mapped) & "pvalue" %in% colnames(deg_mapped)) {
  deg_mapped <- deg_mapped %>%
    mutate(sig = case_when(
      logFC > 1 & pvalue < 0.05 ~ "Up",
      logFC < -1 & pvalue < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    is_gwas = ifelse(mgi_symbol %in% gwas_mouse_genes, "GWAS-linked", "Other"))
  
  png("Volcano_All_DEGs_with_GWAS_highlight.png", width = 2400, height = 2000, res = 300)
  ggplot(deg_mapped, aes(x = logFC, y = -log10(pvalue))) +
    geom_point(aes(color = sig), alpha = 0.5, size = 1.5) +
    geom_point(data = subset(deg_mapped, is_gwas == "GWAS-linked"),
               aes(color = sig, shape = is_gwas), size = 3) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    scale_shape_manual(values = c("GWAS-linked" = 17, "Other" = 16)) + # triangle for GWAS
    theme_minimal(base_size = 14) +
    labs(title = "Volcano Plot (All DEGs, GWAS DEGs highlighted)",
         x = "log2 Fold Change",
         y = "-log10(p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme(legend.position = "top")
  dev.off()
  
  cat("Volcano plot (All DEGs, GWAS DEGs highlighted) saved: Volcano_All_DEGs_with_GWAS_highlight.png\n")
}

# ---- 7. Functional Enrichment ----
if (nrow(obesity_related_degs) > 0) {
  
  entrez_ids <- mapIds(
    org.Mm.eg.db,
    keys = obesity_related_degs$mgi_symbol,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  ) %>% na.omit() %>% unique()
  
  # KEGG
  kegg_res <- enrichKEGG(gene = entrez_ids, organism = 'mmu')
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_df <- as.data.frame(kegg_res)
    write.csv(kegg_df, "KEGG_Enrichment.csv", row.names = FALSE)
    
    # Dotplot
    pdf("KEGG_Dotplot.pdf", width = 8, height = 6)
    dotplot(kegg_res, showCategory = 20) + 
      ggtitle("KEGG Pathway Enrichment (Obesity DEGs)")
    dev.off()
    
    # Barplot
    pdf("KEGG_Barplot.pdf", width = 8, height = 6)
    barplot(kegg_res, showCategory = 20) + 
      ggtitle("KEGG Pathway Enrichment (Obesity DEGs)")
    dev.off()
  }
  
  # Reactome
  reactome_res <- enrichPathway(gene = entrez_ids, organism = "mouse")
  if (!is.null(reactome_res) && nrow(reactome_res) > 0) {
    reactome_df <- as.data.frame(reactome_res)
    write.csv(reactome_df, "Reactome_Enrichment.csv", row.names = FALSE)
    
    # Dotplot
    pdf("Reactome_Dotplot.pdf", width = 8, height = 6)
    dotplot(reactome_res, showCategory = 20) + 
      ggtitle("Reactome Pathway Enrichment (Obesity DEGs)")
    dev.off()
    
    # Barplot
    pdf("Reactome_Barplot.pdf", width = 8, height = 6)
    barplot(reactome_res, showCategory = 20) + 
      ggtitle("Reactome Pathway Enrichment (Obesity DEGs)")
    dev.off()
  }
}
