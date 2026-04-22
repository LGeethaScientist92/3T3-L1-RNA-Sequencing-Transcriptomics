# IL-13 rs20541 RNA-seq Analysis (3T3-L1 Adipocytes)
This repository contains the computational pipeline used for transcriptomic analysis of IL-13 R130 and Q130 variant-treated 3T3-L1 adipocytes.
---
## 📊 Workflow
1. Quality Control – FastQC
2. Transcript Quantification – Salmon (v1.10.2)
3. Gene-level summarisation – tximport
4. Differential Expression – DESeq2
5. Functional Analysis – GO & KEGG enrichment
6. Visualisation – PCA, volcano plots, heatmaps
---
## 📁 Files
### Core pipeline
- `01_fastqc_trim_quant.sh` → Raw FASTQ → Salmon quantification
- `02_deseq2_analysis.R` → DESeq2 analysis + DEG output
- `03_plots_and_enrichment.R` → PCA, heatmaps, volcano, KEGG/GO plots
### Metadata
- `samples.csv` → Sample annotation
---
## Other Scripts:
- `setup_rnaseq_wsl.sh` → Run setup shell
- `TPM Check Script.txt` → QC after running Salmon quantification
- `run_full_pipeline_safe.R` → Full script from Salmon quant.sf files to DEG, GO & KEGG Outputs
- `Generate KEGG Plots.txt/Heatmap.txt/PCA plot.R/functional_figure.R` → Create plots
- `filter_obesity_genes.R` → using Disgenet (steps same as above)
- `Shared pathways.txt` → Plot the shared pathways between Mouse & Human
- `filter_gwas_pipeline.R` → Ensembl → Mouse annotation → GWAS obesity genes → Ortholog mapping → DEG overlap → Obesity-focused enrichment & plots
  ENSEMBL IDs → gene symbols extracted from Mouse Genome Informatics (MGI: https://www.biotools.fr/mouse/ensembl_symbol_converter).
  Mouse phenotype → MGI database: http://www.informatics.jax.org/downloads/reports/ ,filtered for obesity-related Mammalian Phenotype (MP) terms. 
  Obesity-associated human genes → NHGRI-EBI-GWAS Catalog (all associations version 1.0, https://www.ebi.ac.uk/gwas/docs/file-downloads).
  
## 🧬 Experimental Design
| Group     | Samples |
|----------|--------|
| Control  | UnT1, UnT2 |
| R130     | R501, R502 |
| Q130     | Q501, Q502 |
---
## ⚙️ Key Parameters
- Genome: Mus musculus (GRCm39, Ensembl release 111)
- Quantification: Salmon (quasi-mapping mode)
- DEG threshold:
  - adjusted p-value < 0.05
  - |log2FC| > 1
--
## 📈 Outputs
- DEG tables (CSV)
- PCA plots
- Volcano plots
- Top 30 DEG heatmaps
- KEGG / GO enrichment plots
---
## 🌍 Data Availability
RNA-seq data available at:
NCBI Sequence Read Archive (SRA) submission: SUB16136632
---
## 📦 Notes
Additional exploratory scripts 

---
## 👩‍🔬 Author
Geetha Letchumanan  
PhD (Medical Science), Sunway University
