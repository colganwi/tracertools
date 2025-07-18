#!/bin/bash
#SBATCH --job-name=bowtie         # friendly name for job
#SBATCH --nodes=1                 # ensure cpus are on one node
#SBATCH --ntasks=1                # run a single task
#SBATCH --cpus-per-task=4         # number of cpus/threads requested
#SBATCH --mem=4gb                 # memory requested
#SBATCH --partition=20            # partition (queue) to use
#SBATCH --output log/%x_%j.out    # name of output file. %j is jobid
#SBATCH --array=1-2

source activate tracertools

# Parameters
forward_primer="GAATCCAGCTA"
reverse_primer="AGCGGCTAAGG"
r1_pattern="{id}_S{id}_L001_R1_001.fastq.gz"
r2_pattern="{id}_S{id}_L001_R2_001.fastq.gz"
ref="/lab/solexa_weissman/PEtracing_shared/Reference/bowtie/PETS_MFBC_v2_whitelist"
samples_file="./samples.csv"
sites="{"RNF2":332,"HEK3":380,"EMX1":448}"

# Get line from samples file
mapfile -t lines < <(tail -n +2 "$samples_file")
idx=$((SLURM_ARRAY_TASK_ID - 1))
line="${lines[$idx]}"

# parse id and sample
id="$(echo "$line" | cut -d',' -f1 | tr -d '\r\n')"
sample="$(echo "$line" | cut -d',' -f2 | tr -d '\r\n')"
echo "Processing sample: $sample with id: $id"

# format FASTQ files
r1_file="${r1_pattern//\{id\}/$id}"
r2_file="${r2_pattern//\{id\}/$id}"

# Step 1: cutadapt for paired-end reads
mkdir -p ./bam
cutadapt -g "$forward_primer" -G "$reverse_primer" \
  -o "./bam/${id}_R1.trim.fastq" -p "./bam/${id}_R2.trim.fastq" \
  "./fastq/${r1_file}" "./fastq/$r2_file"

# Step 2: bowtie2 alignment for paired-end reads
bowtie2 -x "$ref" \
  -1 "./bam/${id}_R1.trim.fastq" -2 "./bam/${id}_R2.trim.fastq" \
  -S "./bam/${id}.sam" --threads 4 -k 1
rm "./bam/${id}_R1.trim.fastq" "./bam/${id}_R2.trim.fastq"
rm "./bam/${id}.sam"

# Step 3: samtools view, sort, index
samtools view -bS "./bam/${id}.sam" | samtools sort -o "./bam/${sample}.bam"
samtools index "./bam/${sample}.bam"

# Step 4: call alleles
mkdir -p ./data
tracertools alleles-from-bam \
  --bam "./bam/${sample}.bam" \
  --out "./data/${sample}_allele_counts.csv" \
  --barcode_start 270 \
  --barcode_end 300 \
  --site_positions "${sites}" \
  --min_reads 2
