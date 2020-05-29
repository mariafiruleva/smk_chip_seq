rule get_coverage:
    input: bam=rules.sort_bam.output, index_bam=rules.index_bam.output
    output: "coverage/{sample}_{genome}.bw"
    benchmark: "benchmarks/get_coverage/{sample}_{genome}.txt"
    log: "logs/get_coverage/{sample}_{genome}.log"
    conda: "../envs/deeptools.yaml"
    shell: "bamCoverage -b {input.bam} -o {output} &> {log}"