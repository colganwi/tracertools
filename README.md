# TracerTools - a toolkit for for processing PETracer data

## Setup

First time:

```bash
conda env create --file environment.yml
conda activate tracertools
ipython kernel install --user --name tracertools
```

Update:

```bash
conda env update --file environment.yml
```

## Processing 10X data

1. Place .fastq files in a directory named `fastq`.
2. Name .fastq files with sample and library type:
    * `sample1_GEX_S1_R1_001.fastq.gz`
    * `sample1_GEX_S1_R2_001.fastq.gz`
    * `sample1_TS_S1_R1_001.fastq.gz`
    * `sample1_TS_S1_R2_001.fastq.gz`
3. Copy `/templates/cellranger.slurm` to your working directory.
4. Edit the `cellranger.slurm` file to set the sample names.
5. Submit the job array with:
```bash
sbatch cellranger.slurm
```
6. After all jobs are complete, you will have a directory for each sample with the following structure:
```
experiment/
├── bam/
    ├── sample1.bam
    ├── sample1.bam.bai
├── data/
    ├── sample1
        ├── sample1_allele_counts.csv
        ├── sample1_filtered_counts.h5
        ├── sample1_raw_counts.h5
        ├── sample1_summary.html
├── fastq/
    ├── sample1_GEX_S1_R1_001.fastq.gz
    ├── sample1_GEX_S1_R2_001.fastq.gz
    ├── sample1_TS_S1_R1_001.fastq.gz
    ├── sample1_TS_S1_R2_001.fastq.gz
```
7. Make a copy of `/templates/process_10.ipynb` and follow instruction to complete processing and QC.