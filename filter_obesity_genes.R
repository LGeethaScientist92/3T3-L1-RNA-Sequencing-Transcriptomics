library(dplyr)
library(readr)
library(stringr)

# === STEP 1: Load ortholog table (your file) ===
# Ensure file has columns: MGI.symbol, HGNC.symbol
orthologs <- read_csv("human_mouse_orthologs.csv")
cat("Loaded", nrow(orthologs), "ortholog pairs\n")

# === STEP 2: Load local DisGeNET/GWAS reference ===
# Replace with your downloaded file (TSV/CSV). Must have "gene_symbol" and "disease_name" columns.
disgenet <- read_delim("obesity_gwas_genes.csv") 

# Define obesity-related keywords
obesity_keywords <- c("obes", "obesity", "overweight", "bmi", "body mass", "body-mass", "body mass index",
  "adipos", "adipose", "adiposity", "fat mass", "fatness", "adipokine", "leptin",
  "weight gain", "weight loss", "bodyweight", "waist", "waist-hip", "waist to hip", "waist circumference",
  "body fat", "subcutaneous", "visceral",
  "lipid", "lipolysis", "lipogenesis", "triglycerid", "cholesterol", "hdl", "ldl",
  "glucose", "hyperglyc", "hyperglycemia", "insulin resistance", "insulin", "fasting glucose",
  "hba1c", "diabetes", "type 2 diabetes", "t2d", "metabolic syndrome", "metabolic",
  "energy balance", "energy expenditure", "thermogen", "brown adipose", "ucp1",
  "feeding", "hyperphagia", "satiety", "appetite", "ghrelin",
  "metabolism", "body composition", "lean mass", "lipotoxic")

pattern <- str_c(obesity_keywords, collapse = "|")
#pattern <- paste0("(", paste0(keywords, collapse = "|"), ")")
#pattern <- regex(pattern, ignore_case = TRUE)

# Filter DisGeNET/GWAS for obesity associations
obesity_human <- disgenet %>%
  filter(str_detect(tolower(disease_name), pattern)) %>%
  distinct(gene_symbol, disease_name)

cat("Found", nrow(obesity_human), "human obesity-related gene-disease associations\n")

# === STEP 3: Map back to mouse orthologs ===
obesity_orthologs <- obesity_human %>%
  left_join(orthologs, by = c("gene_symbol" = "HGNC.symbol")) %>%
  filter(!is.na(MGI.symbol)) %>%
  distinct()

cat("Mapped to", length(unique(obesity_orthologs$MGI.symbol)), "mouse orthologs\n")


# === STEP 4: Load local MGI phenotype annotations ===
# Download from: http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt
mgi_pheno <- read_delim("MGI_PhenoGenoMP.rpt.txt", col_names = FALSE)
colnames(mgi_pheno) <- c("Allele", "Allele_Symbol", "Background", "MP_ID", "PMID", "MGI_ID") 
mp_terms <- read_delim("VOC_MammalianPhenotype.rpt.txt", delim = "\t", col_names = c("MP_ID", "MP_term", "Definition"))

# Join to add phenotype descriptions
mgi_pheno <- mgi_pheno %>%
  left_join(mp_terms %>% select(MP_ID, MP_term), by = "MP_ID")

# Now you can filter by phenotype keywords
pattern <- "obes|adipos|body mass|bmi|waist|weight|fat|lipid|cholesterol|glucose|insulin|triglyceride|leptin|ghrelin"

mgi_obesity <- mgi_pheno %>%
  filter(str_detect(tolower(MP_term), pattern)) %>%
  distinct(MGI_ID, MP_ID, MP_term)

cat("Found", nrow(mgi_obesity), "mouse obesity-related phenotype associations\n")

# === STEP 2.5: HGNC alias rescue ===
check <- HGNChelper::checkGeneSymbols(unique(obesity_gwas_genes$gene_symbol), unmapped.as.na = FALSE)

obesity_gwas_fixed <- obesity_gwas_genes %>%
  left_join(check %>% select(x, Suggested.Symbol, Approved),
            by = c("gene_symbol" = "x")) %>%
  mutate(final_symbol = ifelse(!is.na(Suggested.Symbol), Suggested.Symbol, gene_symbol))

rescued <- obesity_gwas_fixed %>%
  filter(Approved == FALSE & final_symbol != gene_symbol) %>%
  distinct(gene_symbol, final_symbol)

cat("✅ Rescued", nrow(rescued), "symbols via HGNC alias table\n")

write.csv(rescued, "rescued_gene_symbols.csv", row.names = FALSE)

# Drop junk symbols (intergenic, NR, pseudogenes, etc.)
obesity_gwas_clean <- obesity_gwas_fixed %>%
  filter(!final_symbol %in% c("intergenic", "NR", "NONE")) %>%
  filter(!str_detect(final_symbol, "^MIR")) %>%
  filter(!str_detect(final_symbol, "^LINC")) %>%
  filter(!str_detect(final_symbol, "^LOC"))


# === STEP 5: Report overlap (human + mouse evidence) ===

# Collapse orthologs early → avoid many-to-many join warnings
human_mouse_orthologs_collapsed <- human_mouse_orthologs %>%
  group_by(HGNC_symbol) %>%
  summarise(MGI_symbol = paste(unique(MGI_symbol), collapse = "; "),
            .groups = "drop")

# Map GWAS genes (with rescued names) to mouse orthologs
mapped <- obesity_gwas_clean %>%
  select(gene_symbol, final_symbol, disease_name) %>%
  left_join(human_mouse_orthologs_collapsed,
            by = c("final_symbol" = "HGNC_symbol"))

# Make sure MGI symbols are extracted cleanly from phenotypes
mgi_obesity <- mgi_pheno %>%
  filter(str_detect(tolower(MP_term), pattern)) %>%
  mutate(MGI_symbol = str_extract(Allele_Symbol, "^[^<]+")) %>%  # strip allele part
  distinct(MGI_symbol, MP_ID, MP_term)

# Overlap join (use final_symbol instead of HGNC_symbol)
overlap <- mapped %>%
  inner_join(mgi_obesity, by = "MGI_symbol") %>%
  distinct(final_symbol, MGI_symbol, gene_symbol, disease_name, MP_ID, MP_term)

cat("Overlap genes with both human & mouse obesity evidence:",
    length(unique(overlap$MGI_symbol)), "\n")

# Collapse to avoid row explosions
overlap_collapsed <- overlap %>%
  group_by(final_symbol, MGI_symbol) %>%
  summarise(
    human_symbols = paste(unique(gene_symbol), collapse = "; "),
    disease_names = paste(unique(disease_name), collapse = "; "),
    mp_ids        = paste(unique(MP_ID), collapse = "; "),
    mp_terms      = paste(unique(MP_term), collapse = "; "),
    .groups = "drop"
  )

cat("Collapsed overlap genes:", nrow(overlap_collapsed), "\n")

# === STEP 6: Save outputs ===
write_csv(obesity_orthologs, "obesity_human_mouse_orthologs.csv")
write_csv(mgi_obesity, "obesity_mouse_only.csv")
write_csv(overlap, "obesity_overlap_human_mouse.csv")
write_csv(overlap_collapsed, "obesity_overlap_human_mouse_collapsed.csv")

cat("Results saved:\n - obesity_human_mouse_orthologs.csv\n - obesity_mouse_only.csv\n",
    " - obesity_overlap_human_mouse.csv\n - obesity_overlap_human_mouse_collapsed.csv\n")
