#!/bin/bash

# Step 1: Download datasets
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

# Step 2: Run FastQC
fastqc ACBarrie_R1.fastq.gz ACBarrie_R2.fastq.gz Alsen_R1.fastq.gz Alsen_R2.fastq.gz -o fastqc_results

# Step 3: Trim reads with FastP
fastp -i ACBarrie_R1.fastq.gz -I ACBarrie_R2.fastq.gz -o trimmed_reads/ACBarrie_R1_trimmed.fastq.gz -O trimmed_reads/ACBarrie_R2_trimmed.fastq.gz
fastp -i Alsen_R1.fastq.gz -I Alsen_R2.fastq.gz -o trimmed_reads/Alsen_R1_trimmed.fastq.gz -O trimmed_reads/Alsen_R2_trimmed.fastq.gz

# Step 4: Genome mapping with BWA
bwa index reference.fasta
bwa mem reference.fasta trimmed_reads/ACBarrie_R1_trimmed.fastq.gz trimmed_reads/ACBarrie_R2_trimmed.fastq.gz > ACBarrie_aln.sam
bwa mem reference.fasta trimmed_reads/Alsen_R1_trimmed.fastq.gz trimmed_reads/Alsen_R2_trimmed.fastq.gz > Alsen_aln.sam

# Step 5: Convert SAM to BAM and sort
samtools view -S -b ACBarrie_aln.sam > ACBarrie_aln.bam
samtools sort ACBarrie_aln.bam -o ACBarrie_sorted.bam
samtools view -S -b Alsen_aln.sam > Alsen_aln.bam
samtools sort Alsen_aln.bam -o Alsen_sorted.bam

# Step 6: Variant calling with bcftools
bcftools mpileup -f reference.fasta ACBarrie_sorted.bam | bcftools call -mv -Oz -o ACBarrie_variants.vcf.gz
bcftools mpileup -f reference.fasta Alsen_sorted.bam | bcftools call -mv -Oz -o Alsen_variants.vcf.gz
