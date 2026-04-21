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
## 🧬 Experimental Design
| Group     | Samples |
|----------|--------|
| Control  | UnT1, UnT2 |
| R130     | R501, R502 |
| Q130     | Q501, Q502 |
---
## ⚙️ Key Parameters
- Genome: Mus musculus (GRCm39, Ensembl release 111) Download: https://www.ensembl.org
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
NCBI SRA BioProject: PRJNA1128247
---
## 📦 Notes
Additional exploratory scripts (GWAS integration, ortholog mapping, etc.) are provided in the `/archive` folder.
---
## 👩‍🔬 Author
Geetha Letchumanan  
PhD (Medical Science), Sunway University
