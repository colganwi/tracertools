#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=align          # friendly name for job
#SBATCH --nodes=1                 # ensure cpus are on one node
#SBATCH --ntasks=1                # run a single task
#SBATCH --cpus-per-task=4         # number of cpus/threads requested
#SBATCH --mem=4gb                 # memory requested
#SBATCH --partition=20            # partition (queue) to use
#SBATCH --output log/%x_%j.out    # name of output file. %j is jobid
#SBATCH --array=1-20

# load samples file
samples_file="$1"
if [[ -z "$samples_file" ]]; then
  echo "Usage: $0 <samples.csv>"
  echo "  samples.csv should have lines: prefix,sample"
  exit 1
fi

# specify primers and reference
forward_primer="GAATCCAGCTA"
reverse_primer="AGCGGCTAAGG"
ref="/lab/solexa_weissman/PEtracing_shared/Reference/bowtie/PETS_MFBC_v2_whitelist"

# read samples file into an array (skip header if present)
mapfile -t lines < <(tail -n +2 "$samples_file")

# pick the line corresponding to this array task
idx=$((SLURM_ARRAY_TASK_ID - 1))
line="${lines[$idx]}"

# parse prefix and sample
prefix="$(echo "$line" | cut -d',' -f1)"
sample="$(echo "$line" | cut -d',' -f2)"
echo "Processing sample: $sample with prefix: $prefix"

# locate FASTQ files
r1_file="./fastq/${prefix}_R1.fastq.gz"
r2_file="./fastq/${prefix}_R2.fastq.gz"

# Step 1: cutadapt for paired-end reads
mkdir -p ./bam
cutadapt -g "$forward_primer" -G "$reverse_primer" \
  -o "./bam/${prefix}_R1.trim.fastq" -p "./bam/${prefix}_R2.trim.fastq" \
  "$r1_file" "$r2_file"

# Step 2: bowtie2 alignment
bowtie2 -x "$ref" \
  -1 "./bam/${prefix}_R1.trim.fastq" -2 "./bam/${prefix}_R2.trim.fastq" \
  -S "./bam/${prefix}.sam" --threads 4 -k 1

# clean up trimmed files
rm "./bam/${prefix}_R1.trim.fastq" "./bam/${prefix}_R2.trim.fastq"

# Step 3: samtools view, sort, index
samtools view -bS "./bam/${prefix}.sam" | samtools sort -o "./bam/${sample}.bam"
samtools index "./bam/${sample}.bam"

# Step 4: call alleles
mkdir -p ./data
tracertools alleles-from-bam \
  --bam "./bam/${sample}.bam" \
  --out "./data/${sample}/${sample}_allele_counts.csv" \
  --barcode_start 270 \
  --barcode_end 300 \
  --site_positions "{'RNF2':332,'HEK3':380,'EMX1':448}" \
  --min_reads 2

# final clean up
rm "./bam/${prefix}.sam"