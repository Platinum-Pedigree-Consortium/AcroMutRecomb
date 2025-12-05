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
        "variants/wfmash_raw_var_stats.tsv"
        # expand('variants/{sample}/{sample}.wfmash.sv.vcf', sample=manifest.index),
        # expand('variants/{sample}/{sample}.wfmash.tsv', sample=manifest.index),



rule wgatools_call:
    input:
        paf=find_paf,
        qry= get_query,
        ref=get_target
    output:
        vcf = 'variants/{sample}/{sample}.wfmash.vcf'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "wgatools/1.0.0.1"
    shell:
        """
        samtools faidx {input.qry}
        samtools faidx {input.ref}
        wgatools call -n {wildcards.sample} -f paf --target {input.ref} --query {input.qry} -s -r -o {output.vcf} {input.paf} 
        """

rule separate_sv_snp:
    input:
        vcf = rules.wgatools_call.output.vcf
    output:
        sv = 'variants/{sample}/{sample}.wfmash.sv.vcf',
        snps = 'variants/{sample}/{sample}.wfmash.var.bed'
    params:
        ref_size = get_target_size,
        ref_tig = get_target_tig
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0"
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%END\t[%INFO/SVTYPE]\n" {input.vcf} | awk -v OFS='\\t' '{{if($4=="."){{print $1,$2,$2+1,"SNP"}}else{{print $1,$2,$3,$4}}}}' > {output.snps}
        bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'||SVTYPE=='INV'" -O v {input.vcf} | bcftools filter -e 'REF ~ "N" || ALT ~ "N"' | python /net/eichler/vol28/projects/medical_reference/nobackups/Scripts/Acro/Snakefiles/wgatools_reheader.py {params.ref_tig} {params.ref_size} > {output.sv}
        """

rule vcf2df:
    input:
        vcf = rules.wgatools_call.output.vcf
    output:
        tsv = 'variants/{sample}/{sample}.wfmash.tsv'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        var_list = []
        var_tracker = set()
        for line in open(input.vcf):
            if line.startswith('#'):
                continue
            entries = line.strip().split('\t')

            ref_tig = entries[0]
            ref_sample_name = ref_tig.split('_')[0]
            ref_sample_chr = ref_tig.split('_')[-1]
            ref_sample_hap = ref_tig.split('_')[0] + '_' + ref_tig.split('_')[1]

            qry_tig = entries[9].split(':')[1].split('@')[0]
            qry_pos = int(entries[9].split(':')[1].split('@')[1])
            qry_sample_name = qry_tig.split('_')[0]
            # if qry_sample_name in ['NA12879', 'NA12886']:
            #     continue
            qry_sample_chr = qry_tig.split('_')[-1]
            qry_sample_hap = qry_tig.split('_')[0] + '_' + qry_tig.split('_')[1]


            # if len(entries[3]) >= 50 or len(entries[4]) >= 50:
            #     info_dict = parse_vcf_info_column(entries[7])

            if len(entries[3])==len(entries[4]):
                var_id = f'{qry_sample_name}_{ref_tig}_{entries[1]}_{entries[3]}_{entries[4]}'
                if var_id not in var_tracker:
                    var_list.append([var_id, ref_tig, qry_tig, entries[1], qry_pos, qry_sample_name])
                var_tracker.add(var_id)

        df_vars = pd.DataFrame(var_list, columns=['id', 'ref', 'qry', 'ref_pos', 'qry_pos', 'qry_sample'])
        df_vars.to_csv(output.tsv, sep='\t', header=True, index=False)

def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict

rule raw_var_stats:
    input:
        vcf_list=expand('variants/{sample}/{sample}.wfmash.vcf', sample=manifest.index)
    output:
        stats="variants/wfmash_raw_var_stats.tsv"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        from intervaltree import IntervalTree

        # nucfreq_dir = '/net/eichler/vol28/projects/medical_reference/nobackups/acros/Pedi_vrk22/nucfreq'
        # subregion_dir = '/net/eichler/vol28/projects/medical_reference/nobackups/acros/Pedi_vrk22/pq_contigs_hic/results/pq_haptig_bed'
        #
        # nucfreq_region_dict = {}
        #
        # for line in open('/net/eichler/vol28/projects/medical_reference/nobackups/acros/Pedi_vrk22/pq_contigs_hic/contig_qc/pq_nucfreq_regions.bed'):
        #     tig, start, end, _ = line.strip().split('\t')
        #     if tig not in nucfreq_region_dict:
        #         nucfreq_region_dict[tig] = IntervalTree()
        #     nucfreq_region_dict[tig][int(start): int(end)] = (int(start), int(end))
        #
        # flagger_region_dict = {}
        # for line in open('/net/eichler/vol28/projects/medical_reference/nobackups/acros/Pedi_vrk22/pq_contigs_hic/contig_qc/flagger/acro_flagger_regions.tsv'):
        #     tig, start, end, _ = line.strip().split('\t')
        #     if tig not in flagger_region_dict:
        #         flagger_region_dict[tig] = IntervalTree()
        #     flagger_region_dict[tig][int(start): int(end)] = (int(start), int(end))

        res_list = []
        for vcf_file in input.vcf_list:
            qry_tig = vcf_file.split('/')[1]
            qry_sample_name = qry_tig.split('_')[0]
            qry_sample_chr = qry_tig.split('_')[-1]
            qry_sample_tig = qry_tig.split('_')[2]
            qry_sample_hap = qry_tig.split('_')[0] + '_' + qry_tig.split('_')[1]

            # if qry_sample_name in ['NA12879', 'NA12886']:
            #     continue

            flagged_vars = 0
            var_counter = {'SNP': 0, 'INS': 0, 'DEL': 0, 'INV': 0}
            for line in open(vcf_file):
                if line.startswith('#'):
                    continue
                entries = line.strip().split('\t')

                ref_tig = entries[0]
                ref_sample_name = ref_tig.split('_')[0]
                ref_sample_chr = ref_tig.split('_')[-1]
                ref_sample_hap = ref_tig.split('_')[0] + '_' + ref_tig.split('_')[1]

                qry_tig = entries[9].split(':')[1].split('@')[0]
                qry_pos = int(entries[9].split(':')[1].split('@')[1])
                qry_sample_name = qry_tig.split('_')[0]
                # if qry_sample_name in ['NA12879', 'NA12886']:
                #     continue
                qry_sample_chr = qry_tig.split('_')[-1]
                qry_sample_hap = qry_tig.split('_')[0] + '_' + qry_tig.split('_')[1]
                print(ref_tig, qry_tig)
                if len(entries[3]) >= 50 or len(entries[4]) >= 50:
                    info_dict = parse_vcf_info_column(entries[7])
                    var_counter[info_dict['SVTYPE']] += 1
                else:
                    var_counter['SNP'] += 1

                # query or target sequence in nucfreq flagged regions
                # qry_nucfreq = False if qry_tig not in nucfreq_region_dict else nucfreq_region_dict[qry_tig].overlaps(qry_pos, qry_pos + 1)
                # ref_nucfreq = False if ref_tig not in nucfreq_region_dict else nucfreq_region_dict[ref_tig].overlaps(int(entries[1]), int(entries[1])+ 1)

                # query or target sequence in flagger flagged regions
                # qry_flagger = False if qry_tig not in flagger_region_dict else flagger_region_dict[qry_tig].overlaps(qry_pos,qry_pos + 1)
                # ref_flagger = False if ref_tig not in flagger_region_dict else flagger_region_dict[ref_tig].overlaps(int(entries[1]),int(entries[1]) + 1)


                # if qry_nucfreq or ref_nucfreq or qry_flagger or ref_flagger:
                #     if len(entries[3]) == 1:
                #         flagged_vars += 1
                #     continue

            # if sum(var_counter.values()) > 0:
            res_list.append([qry_sample_name, qry_sample_hap, qry_tig, ref_sample_name, ref_sample_hap, ref_tig, var_counter['SNP'], flagged_vars, var_counter['INS'], var_counter['DEL'], var_counter['INV']])

        res_df = pd.DataFrame(res_list, columns=['qry_sample', 'qry_hap', 'qry_tig', 'ref_sample', 'ref_hap', 'ref_tig', 'Total_SNP', 'InFlag_SNP', 'INS', 'DEL', 'INV'])
        res_df.to_csv(output.stats, sep='\t', header=True, index=False)