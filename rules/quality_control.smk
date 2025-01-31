rule fast_qc:
    input:lambda wildcards: expand("reads/{file_name}", file_name=SAMPLES_INFO.loc[SAMPLES_INFO['SAMPLE'] == wildcards.sample, 'File'])
    output:
          html="qc/fastqc/{sample}.html",
          zip="qc/fastqc/{sample}_fastqc.zip"
    benchmark: "benchmarks/qc/fastqc/{sample}.txt"
    log: "logs/qc/fastqc/{sample}.log"
    threads: 10
    params: ""
    wrapper: "0.59.1/bio/fastqc"

rule multi_qc:
    input: expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['SAMPLE'])
    output: "qc/multiqc/reads.html", directory("qc/multiqc/reads_data")
    benchmark: "benchmarks/qc/multiqc/reads.txt"
    log: "logs/qc/multiqc/reads.log"
    wrapper:
        "0.57.0/bio/multiqc"

rule alignment_qc:
    input: expand("logs/alignment/{sample}_{genome}.log", sample=SAMPLES_INFO['SAMPLE'], genome=config['genome'])
    output: "qc/multiqc/bams.html", directory("qc/multiqc/bams_data")
    benchmark: "benchmarks/qc/multiqc/alignment.txt"
    log: "logs/qc/multiqc/alignment.log"
    conda: "../envs/multi_qc.yaml"
    params:
        outdir="qc/multiqc",
        filename="bams.html"
    shell: "multiqc {input} -o {params.outdir} -n {params.filename}"