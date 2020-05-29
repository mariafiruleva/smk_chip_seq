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

rule sort_bam:
    input: rules.alignment.output
    output: "bams/{sample}_{genome}.sorted.bam"
    benchmark: "benchmarks/sort_bam/{sample}_{genome}.txt"
    log: "logs/sort_bam/{sample}_{genome}.log"
    conda: "../envs/samtools.yaml"
    shell: "samtools sort {input} -o {output} &> {log}"

rule index_bam:
    input: rules.sort_bam.output
    output: "bams/{sample}_{genome}.sorted.bam.bai"
    benchmark: "benchmarks/index_bam/{sample}_{genome}.txt"
    log: "logs/index_bam/{sample}_{genome}.log"
    conda: "../envs/samtools.yaml"
    shell: "samtools index {input} &> {log}"
