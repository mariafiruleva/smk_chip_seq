import re
import pandas as pd

configfile: 'config.yaml'


SAMPLES_INFO = pd.read_csv('data_table.tsv', sep='\t')
SAMPLES_INFO.applymap(lambda x: x.strip() if isinstance(x, str) else x)
SAMPLES_INFO['SAMPLE'] = [f'{SAMPLES_INFO["GSM"][idx]}_{SAMPLES_INFO["Cell"][idx]}_{SAMPLES_INFO["Target"][idx]}' for
                          idx in range(SAMPLES_INFO.shape[0])]

include: "rules/quality_control.smk"
include: "rules/get_reference.smk"
include: "rules/alignment.smk"
include: "rules/get_coverage.smk"

localrules: archive

rule archive:
    input:
         coverage=expand("coverage/{sample}_{genome}.bw", sample=SAMPLES_INFO['SAMPLE'], genome=config['genome']),
         reads_qc="qc/multiqc/reads.html",
         reads_data=directory("qc/multiqc/reads_data"),
         bam_qc="qc/multiqc/bams.html",
         bam_data=directory("qc/multiqc/bams_data")
    output: "chip_seq_results.tar.gz"
    shell: "tar -czvf {output} {input}"


