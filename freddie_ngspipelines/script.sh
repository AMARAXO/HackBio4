
#!/bin/bash

# Directories
DATA_DIR="/home/frederickalloulive/ngs_task/scripts/data"
RESULTS_DIR="/home/frederickalloulive/ngs_task/scripts/results"

# Function to run the pipeline for a dataset
run_pipeline() {
    local dataset_name=$1
    local r1="${DATA_DIR}/${dataset_name}_R1.fastq.gz"
    local r2="${DATA_DIR}/${dataset_name}_R2.fastq.gz"

    echo "Running FastQC for ${dataset_name}..."
    fastqc "$r1" "$r2"

    echo "Trimming reads..."
    fastp -i "$r1" -I "$r2" -o "${RESULTS_DIR}/${dataset_name}_R1_trimmed.fastq.gz" -O "${RESULTS_DIR}/${dataset_name}_R2_trimmed.fastq.gz"

    echo "Indexing reference genome..."
    bwa index "$DATA_DIR/reference.fasta"

    echo "Mapping reads to reference..."
    bwa mem "$DATA_DIR/reference.fasta" "${RESULTS_DIR}/${dataset_name}_R1_trimmed.fastq.gz" "${RESULTS_DIR}/${dataset_name}_R2_trimmed.fastq.gz" > "${RESULTS_DIR}/mapped_${dataset_name}.sam"

    echo "Sorting mapped reads..."
    samtools sort "${RESULTS_DIR}/mapped_${dataset_name}.sam" -o "${RESULTS_DIR}/mapped_${dataset_name}_sorted.bam"

    echo "Indexing sorted BAM file..."
    samtools index "${RESULTS_DIR}/mapped_${dataset_name}_sorted.bam"

    echo "Calling variants..."
    bcftools call -mv -o "${RESULTS_DIR}/called_variants_${dataset_name}.bcf" "${RESULTS_DIR}/mapped_${dataset_name}_sorted.bam"
}

# Run pipeline for each dataset
for dataset in Baxter Drysdale Chara; do
    run_pipeline "$dataset"
done

echo "Pipeline completed successfully!"
