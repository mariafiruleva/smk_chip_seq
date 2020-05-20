import pandas as pd

SAMPLES_INFO = pd.read_csv('data_table.tsv', sep='\t')
SAMPLES_INFO['SAMPLE'] = [f'{SAMPLES_INFO["GSM"][idx]}_{SAMPLES_INFO["Cell"][idx]}_{SAMPLES_INFO["Target"][idx]}' for
                          idx in range(SAMPLES_INFO.shape[0])]

rule fast_qc:
    input: ""
    output: html="qc/fastqc/{sample}.html",
          archive="qc/fastqc/{sample}_fastqc.zip"
    benchmark: "benchmarks/qc/fastqc/{sample}.txt"
    log: "logs/qc/fastqc/{sample}.log"
    conda:
         "envs/fast_qc.yaml"
    shell: ""

rule multi_qc:
    input: expand(rules.fast_qc.output.html)
    output: "qc/multiqc/reads.html"
    benchmark: "benchmarks/qc/multiqc/reads.txt"
    log: "logs/qc/multiqc/reads.log"
    conda:
         "envs/multi_qc.yaml"
    shell: ""
rule get_reference:
    input: ""
    output: ""
    benchmark: "benchmarks/get_reference/get_reference.txt"
    log: "logs/get_reference/get_reference.log"
    shell: ""
rule alignment:
    input: ""
    output: ""
    benchmark: "benchmarks/alignment/alignment.txt"
    log: "logs/alignment/alignment.log"
    shell: ""
rule alignment_qc:
    input: ""
    output: "qc/multiqc/bams.html"
    benchmark: "benchmarks/qc/multiqc/alignment.txt"
    log: "logs/qc/multiqc/alignment.log"
    shell: ""

rule get_coverage:
    input: ""
    output: ""
    benchmark: "benchmarks/get_coverage/get_coverage.txt"
    log: "logs/get_coverage/get_coverage.log"
    shell: ""

rule archive_results:
    input: ""
    output: ""
    benchmark: "benchmarks/archive_results/archive_results.txt"
    log: "logs/archive_results/archive_results.log"
    shell: ""
