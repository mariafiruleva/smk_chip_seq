# smk_chip_seq

## Description

Pipeline takes an input table with samples information, makes quality control using FastQC and summarizes it using Multqc. Then it aligns reads to a UCSC genome build configured in `config.yaml` file, converts BAM files to BigWig files with coverage information and in the end generates final multiqc report including raw qc and bismark qc data.
Input reads are assumed to be single-end.

## Pipeline execution

```
snakemake --use-conda -j 5
```

You can also set number of cores using `-j` parameter.

## Results
[Archive](https://drive.google.com/file/d/1tQN2cxSeIqcRmkX-cJhJ-xTAYVA02rEn/view?usp=sharing)
