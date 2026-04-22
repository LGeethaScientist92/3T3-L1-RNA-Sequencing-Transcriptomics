#!/usr/bin/env bash
set -e

echo "=== Updating system ==="
sudo apt update -y && sudo apt upgrade -y

echo "=== Installing basic tools ==="
sudo apt install -y wget curl unzip git build-essential software-properties-common dirmngr gpg-agent

echo "=== Installing Java (required for Trimmomatic) ==="
sudo apt install -y default-jre

echo "=== Installing libraries needed for R/Bioconductor ==="
sudo apt install -y \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev \
    zlib1g-dev

echo "=== Adding CRAN repo for latest R (4.4) ==="
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys '51716619E084DAB9'
# Replace 'jammy' with 'focal' if your WSL Ubuntu is 20.04
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

echo "=== Installing R 4.4 ==="
sudo apt update
sudo apt install -y r-base r-base-dev

echo "=== Installing FastQC ==="
sudo apt install -y fastqc

echo "=== Installing Trimmomatic ==="
sudo apt install -y trimmomatic

# === Installing Salmon ===
echo "=== Installing Salmon ==="
cd /opt
SALMON_VERSION=1.10.3
wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
tar -xvzf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
mv Salmon-latest_linux_x86_64 salmon-${SALMON_VERSION}
ln -s /opt/salmon-${SALMON_VERSION}/bin/salmon /usr/local/bin/salmon

echo "=== Installing Bioconductor packages in R ==="
Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")'
Rscript -e 'BiocManager::install(version="3.19", ask=FALSE)'
Rscript -e 'BiocManager::install(c("DESeq2", "tximport", "GenomicFeatures", "org.Mm.eg.db", "EnhancedVolcano", "pheatmap"), ask=FALSE)'

echo "=== Setup complete! ==="
echo "✅ R version:" $(R --version | head -n 1)
echo "✅ Salmon version:" $(salmon --version)
echo "✅ FastQC version:" $(fastqc --version)
echo "✅ Trimmomatic version:" $(trimmomatic -version)
