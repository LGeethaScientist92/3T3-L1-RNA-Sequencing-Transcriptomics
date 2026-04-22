#!/bin/bash
set -euo pipefail

# ==========================================
# RNA-seq preprocessing pipeline for 3T3-L1
# FASTQ -> FastQC -> trimming -> Salmon quant
# ==========================================

PROJECT_DIR="/mnt/c/RNASeqProject"
RAW_DIR="$PROJECT_DIR/raw_fastq"
FASTQC_DIR="$PROJECT_DIR/fastqc"
TRIM_DIR="$PROJECT_DIR/trimmed_fastq"
SALMON_INDEX="$PROJECT_DIR/salmon_index"
SALMON_QUANT="$PROJECT_DIR/salmon_quant"
REF_DIR="$PROJECT_DIR/reference"

THREADS=4
TRANSCRIPTOME="$REF_DIR/Mus_musculus.GRCm39.cdna.all.fa"
ADAPTERS="/usr/share/trimmomatic/adapters/TruSeq3-PE.fa"

mkdir -p "$RAW_DIR" "$FASTQC_DIR" "$TRIM_DIR" "$SALMON_INDEX" "$SALMON_QUANT" "$REF_DIR"

echo "=== RNA-seq preprocessing pipeline started ==="

# ==========================================
# Step 0: Check input FASTQ files
# ==========================================
FASTQ_COUNT=$(find "$RAW_DIR" -maxdepth 1 -name "*.fastq.gz" | wc -l)
if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "❌ No FASTQ files found in $RAW_DIR"
    exit 1
fi
echo "✅ Found $FASTQ_COUNT FASTQ files in $RAW_DIR"

# ==========================================
# Step 1: FastQC
# ==========================================
if [ -z "$(find "$FASTQC_DIR" -maxdepth 1 -name "*_fastqc.zip" -print -quit)" ]; then
    echo "=== Running FastQC ==="
    fastqc -t "$THREADS" -o "$FASTQC_DIR" "$RAW_DIR"/*.fastq.gz
else
    echo "✅ FastQC results already exist. Skipping."
fi

# ==========================================
# Step 2: Check adapter file
# ==========================================
if [ ! -f "$ADAPTERS" ]; then
    echo "❌ Adapter file not found: $ADAPTERS"
    echo "Please install Trimmomatic adapters or update the ADAPTERS path."
    exit 1
fi
echo "✅ Adapter file found: $ADAPTERS"

# ==========================================
# Step 3: Check transcriptome FASTA
# ==========================================
if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "❌ Transcriptome FASTA not found: $TRANSCRIPTOME"
    echo "Please download the Mus musculus GRCm39 cDNA FASTA and place it in $REF_DIR"
    exit 1
fi
echo "✅ Transcriptome FASTA found: $TRANSCRIPTOME"

# ==========================================
# Step 4: Trimming with Trimmomatic
# ==========================================
if [ -z "$(find "$TRIM_DIR" -maxdepth 1 -name "*_R1_paired.fq.gz" -print -quit)" ]; then
    echo "=== Running Trimmomatic ==="
    for fq1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
        base=$(basename "$fq1" _R1_001.fastq.gz)
        fq2="$RAW_DIR/${base}_R2_001.fastq.gz"

        if [ ! -f "$fq2" ]; then
            echo "⚠️ Missing pair for $fq1. Skipping."
            continue
        fi

        echo "Processing sample: $base"

        trimmomatic PE -threads "$THREADS" \
            "$fq1" "$fq2" \
            "$TRIM_DIR/${base}_R1_paired.fq.gz" "$TRIM_DIR/${base}_R1_unpaired.fq.gz" \
            "$TRIM_DIR/${base}_R2_paired.fq.gz" "$TRIM_DIR/${base}_R2_unpaired.fq.gz" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
else
    echo "✅ Trimmed FASTQ files already exist. Skipping."
fi

# ==========================================
# Step 5: Build Salmon index
# ==========================================
if [ ! -f "$SALMON_INDEX/versionInfo.json" ]; then
    echo "=== Building Salmon index ==="
    salmon index -t "$TRANSCRIPTOME" -i "$SALMON_INDEX"
else
    echo "✅ Salmon index already exists. Skipping."
fi

# ==========================================
# Step 6: Salmon quantification
# ==========================================
if [ -z "$(find "$SALMON_QUANT" -mindepth 2 -name "quant.sf" -print -quit)" ]; then
    echo "=== Running Salmon quantification ==="
    for fq1 in "$TRIM_DIR"/*_R1_paired.fq.gz; do
        base=$(basename "$fq1" _R1_paired.fq.gz)
        fq2="$TRIM_DIR/${base}_R2_paired.fq.gz"

        if [ ! -f "$fq2" ]; then
            echo "⚠️ Missing trimmed pair for $base. Skipping."
            continue
        fi

        echo "Quantifying sample: $base"

        salmon quant -i "$SALMON_INDEX" -l A \
            -1 "$fq1" -2 "$fq2" \
            -p "$THREADS" \
            -o "$SALMON_QUANT/$base"
    done
else
    echo "✅ Salmon quantification results already exist. Skipping."
fi

echo "🎉 Preprocessing pipeline finished successfully."
echo "Next step: run 02_deseq2_analysis.R"