# tracertools - a toolkit for for processing PETracer data

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

Required software:

- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [samtools](http://www.htslib.org/)
- [slurm](https://slurm.schedmd.com/documentation.html)

## Processing bulk data

1. Place .fastq files in a directory named `fastq`.

```
experiment/
├── fastq/
    ├── sample1_S1_L001_R1.fastq.gz
    ├── sample1_S1_L001_R2.fastq.gz
```

2. Name create CSV sample sheet with id and sample columns:

```
id,sample
1,sample1
```

3. Copy `/templates/bowtie_paired.slurm` or `/templates/bowtie.slurm` to experiment directory.
4. Edit the `.slurm` file to set parameters:
   - `array=1-2` specifies the number of samples to process is parallel (`--array=1-N` for N samples).
   - `r1_pattern` specifies the naming pattern for R1 reads.
   - `r2_pattern` specifies the naming pattern for R2 reads.
   - `forward_primer` and `reverse_primer` are the primers used in the experiment.
   - `ref` is the path to bowtie reference index.
   - `samples_file` is the path to the CSV file with sample information.
   - `sites` is a dictionary with site names and their positions in the reference (set to `"{}"` for intID detection only).
5. Submit the job array with:

```bash
sbatch bowtie_paired.slurm
```

6. After all jobs are complete, you will have a directory for each sample with the following structure:

```
experiment/
├── bam/
    ├── sample1.bam
    ├── sample1.bam.bai
├── data/
    ├── sample1_allele_counts.csv
├── fastq/
    ├── sample1_S1_L001_R1.fastq.gz
    ├── sample1_S1_L001_R2.fastq.gz
```

## Analyzing bulk data

### Number of integrations per sample

```bash
conda activate tracertools
tracertools count-integrations --path ./data/
```

### Edit fraction per sample

```bash
conda activate tracertools
tracertools edit-fractions --path ./data/
```

## Processing 10X data

1. Place .fastq files in a directory named `fastq`.
2. Name .fastq files with sample and library type:

```
experiment/
├── fastq/
    ├── sample1_GEX_S1_R1_001.fastq.gz
    ├── sample1_GEX_S1_R2_001.fastq.gz
    ├── sample1_TS_S1_R1_001.fastq.gz
    ├── sample1_TS_S1_R2_001.fastq.gz
```

3. Copy `/templates/cellranger.slurm` to the experiment directory.
4. Edit the `cellranger.slurm` file to set parameters:
   - `array=1-4` specifies the number of samples to process in parallel (`--array=1-N` for N samples).
   - `samples` is a list of sample names to process.
   - `transcriptome` is the path to the cellranger reference transcriptome.
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
