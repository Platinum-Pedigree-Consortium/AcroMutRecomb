import os.path

import pandas as pd
import numpy as np
import re


configfile: "liftover.yaml"


MANIFEST = config.get("MANIFEST")
BED = config.get("BED")
bed_df = pd.read_csv(
    BED, sep="\t", header=None, names=["chr", "start", "end"], dtype=str
)
REF = config['REF']
bed_df["NAME"] = bed_df["chr"] + "_" + bed_df["start"] + "_" + bed_df["end"]

bed_df.set_index("NAME", inplace=True)

manifest_df = pd.read_csv(MANIFEST, sep='\t', index_col=['SAMPLE'])



# rustybam/0.1.16

def cigar_tuple(cigar):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]
    tuples = []
    for i in range(len(lengths)):
        tuples.append([int(lengths[i]), ops[i]])
    return tuples

def find_asm(wildcards):
    return manifest_df.at[wildcards.sample, "HAP"]

def find_paf(wildcards):
    return manifest_df.at[wildcards.sample, "PAF"]

def find_tigs(wildcards):
    tig_dict = {}
    for sample in manifest_df.index:
        for region in bed_df.index:
            tig_dict[f'{sample}_{region}_pq'] = rules.tag_contigs.output.pq.format(sample=sample,region=region)
            tig_dict[f'{sample}_{region}_p'] = rules.tag_contigs.output.p.format(sample=sample,region=region)
    return tig_dict

wildcard_constraints:
    sample="|".join(manifest_df.index),
    region="|".join(bed_df.index),

localrules:
    liftover,
    # moddotplot,


rule liftover:
    input:
        # expand("results/asm_to_chm13/{sample}.paf",sample=manifest_df.index),
        # expand("results/{region}/contigs/{sample}_pq_contig.bed",sample=manifest_df.index,region=bed_df.index),
        # expand("results/{region}/fasta/{sample}_pq_contig.fa",sample=manifest_df.index,region=bed_df.index),
        expand("results/{region}/moddotplot/{sample}/success", sample=manifest_df.index, region=bed_df.index),
        # "results/lifted_contigs.tab",

# rule moddotplot:
#     input:
#         expand("results/{region}/moddotplot/{sample}_pq_contig.png",sample=manifest_df.index,region=bed_df.index),

rule subset_target_region:
    input:
        bed=BED,
    output:
        bed="tmp/{region}.bed",
    resources:
        mem=10,
        hrs=24,
        disk_free=1
    run:
        out_df = bed_df.loc[wildcards.region]
        with open(output.bed, "w") as out_file:
            out_file.write("\t".join(out_df[["chr", "start", "end"]]))
            out_file.write("\n")


rule mm2_aln:
    input:
        query=find_asm,
    output:
        aln='results/asm_to_chm13/{sample}.paf',
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "minimap2/2.28",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        minimap2 -x asm20 --secondary=no -s 25000 -K 8G -c -t {threads} --eqx --cs {REF} {input.query} > {output.aln}
        """


rule liftoff:
    input:
        bed=rules.subset_target_region.output.bed,
        paf='results/asm_to_chm13/{sample}.paf',
    output:
        paf="results/{region}/paf/{sample}.paf",
    resources:
        mem=24,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        mkdir -p $( dirname {output.paf} )
        rustybam liftover --bed {input.bed} {input.paf} > {output.paf}
        """


rule trim_paf:
    input:
        paf="results/{region}/paf/{sample}.paf",
    output:
        paf="results/{region}/paf/{sample}_trimmed.paf",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        rustybam trim-paf {input.paf} > {output.paf}
        rm {input.paf}
        """

rule paf_stats:
    input:
        paf="results/{region}/paf/{sample}_trimmed.paf",
    output:
        stats="results/{region}/paf/{sample}_trimmed.stats",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        rustybam stats --paf --qbed {input.paf} > {output.stats}
        """

rule tag_contigs:
    input:
        paf = "results/{region}/paf/{sample}_trimmed.paf"
    output:
        pq = "results/{region}/contigs/{sample}_pq_contig.bed",
        p = "results/{region}/contigs/{sample}_p_contig.bed",
        aln = "results/{region}/contigs/{sample}_contig_aln.paf"

    threads: 1
    resources:
        mem=10,
        hrs=24,
        disk_free=1
    run:
        centro_dict = {
                'chr13': {
                    'p' : (15547593, 16522942),
                    'q' : (16522942, 17498291),
                },
                'chr14': {
                    'p' : (10092112, 11400261),
                    'q' : (11400261, 12708411),
                },
                'chr15': {
                    'p' : (16678794, 17186630),
                    'q' : (17186630, 17694466),
                },
                'chr21': {
                    'p' : (10962853, 11134529),
                    'q' : (11134529, 11306205),
                },
                'chr22': {
                    'p' : (12788180, 14249622),
                    'q' : (14249622, 15711065),
                }
            }

        seq2len_dict = {}
        query2target2info_dict = {}
        ref_tag = wildcards.region.split('_')[0]
        sample_name = wildcards.sample.split('_')[0]
        qry_chrom_tracker = {}
        with open(input.paf) as f:
            for line in f:
                query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]

                cg = [i.split(":")[-1] for i in line.strip().split("\t")[12:] if i[:2] == 'cg'][0]
                cg_tuple = cigar_tuple(cg)
                iden = round((sum([int(i[0]) for i in cg_tuple if i[1] == '=' or i[1] == 'M']) / sum(
                    [int(i[0]) for i in cg_tuple if i[1] in {'=', 'M', 'X', 'D', 'I'}])) * 100,2)

                query_len = int(query_len)
                query_start = int(query_start)
                query_end = int(query_end)
                target_len = int(target_len)
                target_start = int(target_start)
                target_end = int(target_end)

                seq2len_dict[query] = query_len
                seq2len_dict[target] = target_len

                if query not in qry_chrom_tracker:
                    qry_chrom_tracker[query] = []
                qry_chrom_tracker[query].append(target)


                if int(alignment_len) > 100000 and iden >= 90:
                    if query not in query2target2info_dict:
                        query2target2info_dict[query] = {'p': [], 'q': []}

                    if target_end - 1000000 >= centro_dict[target]['q'][1] and int(alignment_len) > 1000000:
                        arm = 'q'
                        query2target2info_dict[query][arm].append((
                        (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                        line.strip()))

                    if target_start <= centro_dict[target]['p'][1]:
                        arm = 'p'
                        query2target2info_dict[query][arm].append((
                        (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                        line.strip()))


        query2arm2target2info_dict = {}

        fout_pq = open(output.pq,'w')
        fout_p = open(output.p, 'w')
        faln = open(output.aln, 'w')
        pq_contig_tracker = []
        for query, query2arm_dict in query2target2info_dict.items():
            query_len = seq2len_dict[query]
            ## contigs aligned to both p and q arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) > 0:
                distal_aln = sorted(query2arm_dict['q'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                # iden = distal_aln[5]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                target_chrom = distal_aln[3]

                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{query}_{target_chrom}',file=fout_pq)
                ## used for Arang's annotated assembly
                new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{new_query}',file=fout_pq)
                # pq_contig_tracker.append([query, qry_seq_start, qry_seq_end, f'{wildcards.sample}_{query}_{target_chrom}', iden])

            ## contig only aligned to p arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) == 0:
                distal_aln = sorted(query2arm_dict['p'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                target_chrom = distal_aln[3]
                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{query}_{target_chrom}',file=fout_p)
                ## used for Arang's annotated assembly
                new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{new_query}',file=fout_p)

            for arm, qry2arm_list in query2arm_dict.items():
                for aln_info in qry2arm_list:
                    print(aln_info[-1], file=faln)


        # if len(pq_contig_tracker) > 1:
        #     sorted_pq_contig = sorted(pq_contig_tracker, key=lambda x:x[-1], reverse=True)[0]
        #     print(f'{sorted_pq_contig[0]}\t{sorted_pq_contig[1]}\t{sorted_pq_contig[2]}\t{sorted_pq_contig[3]}', file=fout_pq)
        # else:
        #     print(f'{pq_contig_tracker[0][0]}\t{pq_contig_tracker[0][1]}\t{pq_contig_tracker[0][2]}\t{pq_contig_tracker[0][3]}',file=fout_pq)



rule group_tigs:
    input:
        unpack(find_tigs),
    output:
        tigs = "results/lifted_contigs.tsv",
        tab = "results/lifted_contigs.tab",
    threads: 1
    resources:
        mem=10,
        hrs=24,
        disk_free=1

    run:
        acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
        fout = open(output.tigs, 'w')
        print('Sample\tArm\tTig\tTig_start\tTig_end\tSize\tChrom',file=fout)
        haps = open(output.tab, 'w')
        print('SAMPLE\tHAP', file=haps)

        for sample in manifest_df.index:
            qry_chrom_tracker = {}
            qrylen_dict = {}
            with open(f'results/asm_to_chm13/{sample}.paf') as f:
                for line in f:
                    query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]
                    qrylen_dict[query] = query_len
                    if int(alignment_len) < 200000:
                        continue
                    if query not in qry_chrom_tracker:
                        qry_chrom_tracker[query] = []
                    qry_chrom_tracker[query].append(target)

            for region in bed_df.index:
                for line in open(input[f'{sample}_{region}_pq']):
                    entries = line.strip().split('\t')
                    qlen = qrylen_dict[entries[0]]
                    print(f"{sample}\tpq\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)
                    if os.path.exists(f'results/{region}/fasta/{sample}_pq_contig.fa'):
                        print(f'{entries[3]}\tresults/{region}/fasta/{sample}_pq_contig.fa', file=haps)

                for line in open(input[f'{sample}_{region}_p']):
                    entries = line.strip().split('\t')
                    qlen = qrylen_dict[entries[0]]
                    not_acro = False
                    for ele in qry_chrom_tracker[entries[0]]:
                        if ele not in acros:
                            not_acro = True
                            break
                    if not not_acro:
                        print(f"{sample}\tp\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)

# rule get_pq_fa:
#     input:
#         bed="results/{region}/contigs/{sample}_pq_contig.bed",
#         hap = find_asm,
#     output:
#         fa="results/{region}/fasta/{sample}_pq_contig.fa"
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "miniconda/4.12.0",
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1,
#     threads: 1
#     shell:
#         """
#         seqtk subseq {input.hap} {input.bed} > {output.fa}
#         """
rule get_pq_fa:
    input:
        bed="results/{region}/contigs/{sample}_pq_contig.bed",
        hap = find_asm,
    output:
        fa="results/{region}/fasta/{sample}_pq_contig.fa"
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        if [ -s {input.bed} ]; then
            header=$(awk '{{print $4}}' OFS='' FS='\\t' {input.bed})
            seqtk subseq {input.hap} {input.bed} | sed "s/^>.*/>$header/" > {output.fa}
        fi
        """

# bedtools getfasta -fi {input.hap} -bed {input.bed} -fo {output.fa} -nameOnly

rule pq_selfplot:
    input:
        bed = "results/{region}/contigs/{sample}_pq_contig.bed",
        fa = "results/{region}/fasta/{sample}_pq_contig.fa"
    output:
        out = "results/{region}/moddotplot/{sample}/success"

    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
        "moddotplot/0.9.0"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        if [ -s {input.bed} ]; then
            moddotplot static -f {input.fa} --no-hist -o $( dirname {output.out} ) --compare
        fi
        echo Finished! > {output.out} 
        """