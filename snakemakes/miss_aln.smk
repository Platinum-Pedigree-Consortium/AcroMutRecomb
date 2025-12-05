import pandas as pd


manifest = pd.read_csv('detected_var.tab', sep='\t', index_col=['SAMPLE'])

def get_query(wildcards):
    return manifest.at[wildcards.sample, 'QUERY_FA']

def get_target(wildcards):
    return manifest.at[wildcards.sample, 'TARGET_FA']

def find_paf(wildcards):
    return manifest.at[wildcards.sample, 'PAF']

def get_target_tig(wildcards):
    return manifest.at[wildcards.sample, 'TARGET'].replace(';', ',')

def get_target_size(wildcards):
    return manifest.at[wildcards.sample, 'TARGET_SIZE'].replace(';', ',')

rule all:
    input:
        expand('variants/{sample}/slider_align_w10k.failed.wfmash.bed',sample=manifest.index),

rule slider_align:
    input:
        paf=find_paf,
    output:
        stats='variants/{sample}/slider_align_w10k.wfmash.bed'
    params:
        ref = get_target_tig,
        size= get_target_size
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.33"
    shell:
        """
        rb liftover --bed <(bedtools makewindows -w 10000 -b <(printf "{params.ref}\t0\t{params.size}\n")) {input.paf} | rb stats --paf > {output.stats}
        """

rule slider_failed:
    input:
        stats='variants/{sample}/slider_align_w10k.wfmash.bed'
    output:
        bed='variants/{sample}/slider_align_w10k.failed.wfmash.bed'
    params:
        ref = get_target_tig,
        size= get_target_size
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "rustybam/0.1.33"
    shell:
        """
        awk '$10<99.99' {input.stats} | cut -f 1-3 | sortBed -i /dev/stdin | bedtools merge -i /dev/stdin > {output.bed}
        """