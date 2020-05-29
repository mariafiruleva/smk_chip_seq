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
    conda: '../envs/bowtie.yaml'
    threads: 4
    params:
        basename=f'indexes/{config["genome"]}/{config["genome"]}'
    shell: "bowtie2-build {input} {params.basename} &> {log}"
