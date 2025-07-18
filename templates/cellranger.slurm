#!/bin/bash
#SBATCH --job-name=cellranger          # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=32             # number of cpus/threads requested.
#SBATCH --mem=128gb                    # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
#SBATCH --output=cellranger-%j.out     # name of output file.  %j is jobid
#SBATCH --array=1-4                    # array of jobs to submit

source activate tracertools

# Update list of samples
samples=(
  "embryo27_t1"
  "embryo27_t2"
  "embryo28_t1"
  "embryo28_t2"
)

# Directory containing the fastq files
fastq_dir="./fastq"
fastq_dir=$(realpath "$fastq_dir")

# Path to the transcriptome reference
transcriptome="/lab/solexa_weissman/PEtracing_shared/Reference/cellranger/mm10_PETS_v2"
transcriptome=$(realpath "$transcriptome")


# Temarborium for the libraries.csv content
libraries_temarborium="fastqs,sample,library_type
{fastq_dir},{sample}_TS,Gene Expression
{fastq_dir},{sample}_GEX,Gene Expression"

#Make output directories
mkdir -p ./bam
mkdir -p ./data

# Sample
sample=${samples[($SLURM_ARRAY_TASK_ID-1)]}
echo "Processing sample: ${sample}"

# Run cellranger count 
libraries_csv="${sample}_libraries.csv"
echo "Running cellranger count for sample: ${sample}"
cellranger count \
  --id="$sample" \
  --libraries="$libraries_csv" \
  --transcriptome="$transcriptome" \
  --localmem=120 \
  --nosecondary

# Clean up
rm "$libraries_csv"
output_dir="./$sample/outs"
data_dir="./data/${sample}"
mkdir -p $data_dir
# Bam
bam_file="${output_dir}/possorted_genome_bam.bam"
mv "$bam_file" "./bam/${sample}.bam"
bai_file="${output_dir}/possorted_genome_bam.bam.bai"
mv "$bai_file" "./bam/${sample}.bai"
# Counts
filtered_h5="${output_dir}/filtered_feature_bc_matrix.h5"
mv "$filtered_h5" "${data_dir}/${sample}_filtered_counts.h5"
raw_h5="${output_dir}/raw_feature_bc_matrix.h5"
mv "$raw_h5" "${data_dir}/${sample}_raw_counts.h5"
# Summary
summary="${output_dir}/web_summary.html"
mv "$summary" "${data_dir}/${sample}_summary.html"
# Remove output directory
rm -r "./${sample}"

# Call TS alleles
echo "Calling TS alleles for sample: ${sample}"
tracertools alleles-from-bam \
  --bam "./bam/${sample}.bam" \
  --out "./data/"${sample}"/${sample}"_allele_counts.csv \
  --barcode_start 270 \
  --barcode_end 300 \
  --site_positions  "{'RNF2':332,'HEK3':380,'EMX1':448}" \
  --min_reads 2