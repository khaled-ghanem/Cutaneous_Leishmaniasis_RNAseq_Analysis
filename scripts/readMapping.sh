#!/bin/bash

# Set up directories
DATA_DIR="./data"
FASTQ_DIR="$DATA_DIR/fastq"
RESULTS_DIR="./results"

# Create results subdirectories
mkdir -p "$RESULTS_DIR/fastqc"
mkdir -p "$RESULTS_DIR/kallisto"
mkdir -p "$RESULTS_DIR/multiqc"



# Create and activate conda environment
echo "Creating conda environment..."
conda create -n rna_mapping -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    kallisto=0.50.1 \
    multiqc=1.21 \
    -y

echo "Activating conda environment..."
source /home/khaled/miniconda3/etc/profile.d/conda.sh
conda activate rna_mapping

# Check if activation was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment"
    exit 1
fi

# Navigate to fastq directory
pwd
cd "./data/fastq" || exit 1

# Run fastqc to check the quality of our fastq files
echo "Running FastQC..."
fastqc *.gz -t 4 -o "../../results/fastqc"

# Navigate to data directory for reference genome
pwd
cd "../" || exit 1

# Build an index from the reference fasta file
echo "Building Kallisto index..."
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa
cd "./fastq"

# Map reads to the indexed reference host transcriptome
echo "Mapping reads with Kallisto..."

# Healthy subjects (HS)
echo "Processing healthy subjects..."
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/HS01" -t 4 --single -l 250 -s 30 SRR8668755.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/HS02" -t 4 --single -l 250 -s 30 SRR8668756.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/HS03" -t 4 --single -l 250 -s 30 SRR8668757.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/HS04" -t 4 --single -l 250 -s 30 SRR8668758.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/HS05" -t 4 --single -l 250 -s 30 SRR8668759.fastq.gz 

# Cutaneous leishmaniasis (CL) patients
echo "Processing CL patients..."
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/CL08" -t 4 --single -l 250 -s 30 SRR8668769.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/CL10" -t 4 --single -l 250 -s 30 SRR8668771.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/CL11" -t 4 --single -l 250 -s 30 SRR8668772.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/CL12" -t 4 --single -l 250 -s 30 SRR8668773.fastq.gz 
kallisto quant -i "../Homo_sapiens.GRCh38.cdna.all.index" -o "../../results/kallisto/CL13" -t 4 --single -l 250 -s 30 SRR8668774.fastq.gz 

# Summarize fastqc and kallisto mapping results using MultiQC
echo "Running MultiQC..."
pwd
cd "../../results" || exit 1
multiqc . -o "./multiqc"

# Deactivate conda environment
conda deactivate

echo "Finished! Results are in: $RESULTS_DIR"
echo "  - FastQC reports: $RESULTS_DIR/fastqc"
echo "  - Kallisto results: $RESULTS_DIR/kallisto"
echo "  - MultiQC summary: $RESULTS_DIR/multiqc"
