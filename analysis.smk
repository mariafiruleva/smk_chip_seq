import re

configfile: 'config.yaml'

import pandas as pd

SAMPLES_INFO = pd.read_csv('data_table.tsv', sep='\t')
SAMPLES_INFO.applymap(lambda x: x.strip() if isinstance(x, str) else x)
SAMPLES_INFO['SAMPLE'] = [f'{SAMPLES_INFO["GSM"][idx]}_{SAMPLES_INFO["Cell"][idx]}_{SAMPLES_INFO["Target"][idx]}' for
                          idx in range(SAMPLES_INFO.shape[0])]

rule archive:
    input:
         coverage=expand("coverage/{sample}_{genome}.bw", sample=SAMPLES_INFO['SAMPLE'], genome=config['genome']),
         reads_qc="qc/multiqc/reads.html",
         reads_data=directory("qc/multiqc/reads_data"),
         bam_qc="qc/multiqc/bams.html",
         bam_data=directory("qc/multiqc/bams_data")
    output: "chip_seq_results.tar.gz"
    shell: "tar -czvf {output} {input}"

rule fast_qc:
    input:lambda wildcards: expand("reads/{file_name}", file_name=SAMPLES_INFO.loc[SAMPLES_INFO['SAMPLE'] == wildcards.sample, 'File'])
    output:
          html="qc/fastqc/{sample}.html",
          archive="qc/fastqc/{sample}_fastqc.zip"
    benchmark: "benchmarks/qc/fastqc/{sample}.txt"
    log: "logs/qc/fastqc/{sample}.log"
    conda:
         "envs/fast_qc.yaml"
    threads: 10
    params:
          out_dir="qc/fastqc",
          out_sample=lambda wildcards: re.sub('SRX\d*/|.fastq', '', list(
              SAMPLES_INFO.loc[SAMPLES_INFO['SAMPLE'] == wildcards.sample, 'File'].to_dict().values())[0])
    shell: """
    fastqc {input} -o {params.out_dir} 2> {log}
    mv {params.out_dir}/{params.out_sample}*zip {output.archive}
    mv {params.out_dir}/{params.out_sample}*html {output.html}
    """

rule multi_qc:
    input: expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['SAMPLE'])
    output: "qc/multiqc/reads.html", directory("qc/multiqc/reads_data")
    benchmark: "benchmarks/qc/multiqc/reads.txt"
    log: "logs/qc/multiqc/reads.log"
    wrapper:
        "0.57.0/bio/multiqc"
rule get_reference:
    output: "indexes/{genome}/{genome}.fa.gz"
    benchmark: "benchmarks/get_reference/get_reference_{genome}.txt"
    log: "logs/get_reference/get_reference_{genome}.log"
    shell:
         "wget -O {output} http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz &> {log}"

rule bowtie2_index:
    input: rules.get_reference.output
    output:
          'indexes/{genome}/{genome}.1.bt2', 'indexes/{genome}/{genome}.2.bt2',
          'indexes/{genome}/{genome}.3.bt2', 'indexes/{genome}/{genome}.4.bt2',
           'indexes/{genome}/{genome}.rev.1.bt2', 'indexes/{genome}/{genome}.rev.2.bt2'
    log: "logs/build/{genome}.log"
    benchmark: "benchmarks/build/{genome}.txt"
    conda: 'envs/bowtie.yaml'
    threads: 4
    params:
        basename=f'indexes/{config["genome"]}/{config["genome"]}'
    shell: "bowtie2-build {input} {params.basename} &> {log}"

rule alignment:
    input: sample=lambda wildcards: expand("reads/{file_name}", file_name=SAMPLES_INFO.loc[SAMPLES_INFO['SAMPLE'] == wildcards.sample, 'File']),
           index_path=rules.bowtie2_index.output
    output: temp("bams/{sample}_{genome}.bam")
    benchmark: "benchmarks/alignment/{sample}_{genome}.txt"
    log: "logs/alignment/{sample}_{genome}.log"
    params:
        index='indexes/{genome}/{genome}',
        extra=''
    threads: 8
    wrapper:
        "0.57.0/bio/bowtie2/align"
rule alignment_qc:
    input: expand("logs/alignment/{sample}_{genome}.log", sample=SAMPLES_INFO['SAMPLE'], genome=config['genome'])
    output: "qc/multiqc/bams.html", directory("qc/multiqc/bams_data")
    benchmark: "benchmarks/qc/multiqc/alignment.txt"
    log: "logs/qc/multiqc/alignment.log"
    conda: "envs/multi_qc.yaml"
    params:
        outdir="qc/multiqc",
        filename="bams.html"
    shell: "multiqc {input} -o {params.outdir} -n {params.filename}"


rule sort_bam:
    input: rules.alignment.output
    output: "bams/{sample}_{genome}.sorted.bam"
    benchmark: "benchmarks/sort_bam/{sample}_{genome}.txt"
    log: "logs/sort_bam/{sample}_{genome}.log"
    conda: "envs/samtools.yaml"
    shell: "samtools sort {input} -o {output} &> {log}"

rule index_bam:
    input: rules.sort_bam.output
    output: "bams/{sample}_{genome}.sorted.bam.bai"
    benchmark: "benchmarks/index_bam/{sample}_{genome}.txt"
    log: "logs/index_bam/{sample}_{genome}.log"
    conda: "envs/samtools.yaml"
    shell: "samtools index {input} &> {log}"

rule get_coverage:
    input: bam=rules.sort_bam.output, index_bam=rules.index_bam.output
    output: "coverage/{sample}_{genome}.bw"
    benchmark: "benchmarks/get_coverage/{sample}_{genome}.txt"
    log: "logs/get_coverage/{sample}_{genome}.log"
    conda: "envs/deeptools.yaml"
    shell: "bamCoverage -b {input.bam} -o {output} &> {log}"
