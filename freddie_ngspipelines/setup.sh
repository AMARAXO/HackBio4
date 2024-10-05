


#!/bin/bash

echo "Installing required tools..."

# Install dependencies
sudo apt-get update
sudo apt-get install -y fastqc fastp bwa samtools bcftools

echo "All tools installed!"
