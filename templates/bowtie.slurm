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
ref="/lab/solexa_weissman/PEtracing_shared/Reference/bowtie/PETS_MFBC_v2_whitelist"
samples_file="./samples.csv"
sites="{"RNF2":332}"

# Get line from samples file
mapfile -t lines < <(tail -n +2 "$samples_file")
idx=$((SLURM_ARRAY_TASK_ID - 1))
line="${lines[$idx]}"

# parse id and sample
id="$(echo "$line" | cut -d',' -f1 | tr -d '\r\n')"
sample="$(echo "$line" | cut -d',' -f2 | tr -d '\r\n')"
echo "Processing sample: $sample with id: $id"

# Format FASTQ file name
r1_file="${r1_pattern//\{id\}/$id}"

# Step 1: cutadapt
mkdir -p ./bam
cutadapt -g "$forward_primer" -a "$reverse_primer" --action retain "./fastq/${r1_file}" > "./bam/${id}.trim.fastq"

# Step 2: bowtie2 alignment
bowtie2 -x ${ref} -U "./bam/${id}.trim.fastq" -S "./bam/${id}.sam" --threads 4 -k 1
rm "./bam/${id}.trim.fastq"

# Step 3: samtools view, sort, index
samtools view -bS "./bam/${id}.sam" | samtools sort -o "./bam/${sample}.bam"
samtools index "./bam/${sample}.bam"
rm "./bam/${id}.sam"

# Step 4: call alleles
mkdir -p ./data
tracertools alleles-from-bam \
  --bam "./bam/${sample}.bam" \
  --out "./data/${sample}_allele_counts.csv" \
  --barcode_start 270 \
  --barcode_end 300 \
  --site_positions "${sites}" \
  --min_reads 2

