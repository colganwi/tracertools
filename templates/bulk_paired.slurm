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
  echo "  samples.csv should have lines: sample_prefix,sample"
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
sample_prefix="$(echo "$line" | cut -d',' -f1)"
sample="$(echo "$line" | cut -d',' -f2)"

# locate FASTQ files
r1_file="../PETS/${sample_prefix}_L001_R1_001.fastq.gz"
r2_file="../PETS/${sample_prefix}_L001_R2_001.fastq.gz"

# Step 1: cutadapt for paired-end reads
cutadapt -g "$forward_primer" -G "$reverse_primer" \
  -o "${sample_prefix}_R1.trim.fastq" -p "${sample_prefix}_R2.trim.fastq" \
  "$r1_file" "$r2_file"

# Step 2: bowtie2 alignment
bowtie2 -x "$ref" \
  -1 "${sample_prefix}_R1.trim.fastq" -2 "${sample_prefix}_R2.trim.fastq" \
  -S "${sample_prefix}.sam" --threads 4 -k 1

# clean up trimmed files
rm "${sample_prefix}_R1.trim.fastq" "${sample_prefix}_R2.trim.fastq"

# Step 3: samtools view, sort, index
samtools view -bS "${sample_prefix}.sam" \
  | samtools sort -o "${sample}.sorted.bam"
samtools index "${sample}.sorted.bam"

# final clean up
rm "${sample_prefix}.sam"