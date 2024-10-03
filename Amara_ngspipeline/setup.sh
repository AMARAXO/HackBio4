#!/bin/bash
# Install necessary tools for NGS pipeline
sudo apt-get update

# Install wget to download datasets
sudo apt-get install -y wget

# Install FastQC for quality control
sudo apt-get install -y fastqc

# Install FastP for trimming reads
sudo apt-get install -y fastp

# Install BWA for genome mapping
sudo apt-get install -y bwa

# Install samtools and bcftools for variant calling
sudo apt-get install -y samtools bcftools

