#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/4/25
'''
import os
import pandas as pd
from intervaltree import IntervalTree

def get_overlaps(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))
def get_valid_dnms(var_dir, aligner, children, fp_trans):

    classified_tbl = pd.read_csv('./classified_contig_matches.txt',sep='\t',header=[0])
    classified_tbl = classified_tbl.loc[~classified_tbl['G3-TIG'].isin(fp_trans)]

    classified_tbl['G3-SAMPLE'] = classified_tbl.apply(lambda row:row['G3-TIG'].split('_')[0], axis=1)
    contig_size = pd.read_csv('./all_pq_contigs.size', sep='\t',usecols=[0, 1], names=['contig', 'size'])
    contig_size.set_index('contig', inplace=True)
    qarm_range = pd.read_csv('./contig_qc_summary.tsv', sep='\t', header=0)
    qarm_range = qarm_range.loc[qarm_range['HAP']!='.']
    qarm_range['size'] = qarm_range.apply(lambda row: contig_size.at[row['HAP'], 'size'] - int(row['qARM_START']), axis=1)
    qarm_range['asat_size'] = qarm_range.apply(lambda row: contig_size.at[row['HAP'], 'size'] - int(row['aSat_START']), axis=1)

    qarm_range['end'] = qarm_range.apply(lambda row: contig_size.at[row['HAP'], 'size'], axis=1)
    qarm_range.set_index('HAP', inplace=True)
    # qarm_range['size'] = qarm_range.apply(lambda row:row['end'] - row['start'], axis=1)
    # qarm_range.set_index('tig_name', inplace=True)

    flagged_dict = {}
    for line in open(f'./flagger_nucfreq_merged.bed'):
    # for line in open(f'{VOL28}/acros/Pedi_vrk22/refOri_pq_contigs/contig_qc/pq_nucfreq_regions.bed'):
    #     chrom, start, end, _, _ = line.strip().split('\t')
        chrom, start, end = line.strip().split('\t')
        if chrom not in flagged_dict:
            flagged_dict[chrom] = IntervalTree()
        flagged_dict[chrom][int(start): int(end)] = [int(start), int(end)]

    mother_dnms = []
    fater_dnms = []
    mut_rate = []
    mut_chrom = []
    callable_stats = []

    if not os.path.exists(f'{var_dir}/{aligner}_val/validated_calls/region_passed_true_dnm'):
        os.mkdir(f'{var_dir}/{aligner}_val/validated_calls/region_passed_true_dnm')

    for sample in children:
        g3_tig_matches = classified_tbl.loc[classified_tbl['G3-SAMPLE'] == sample]
        de_novo = {'father': [0, 0, 0, 0, 0], 'mother': [0, 0, 0, 0, 0]}

        ## short arm bases include aSat
        bases = {'father': 0, 'mother': 0}
        ## short arm bases exclude aSat
        bases_ex_aSat = {'father': 0, 'mother': 0}

        bases_chrom = {'father': {'chr13':0, 'chr14':0, 'chr15':0, 'chr21':0, 'chr22':0}, 'mother': {'chr13':0, 'chr14':0, 'chr15':0, 'chr21':0, 'chr22':0}}

        # for i, chrom in enumerate(['chr13', 'chr14', 'chr15', 'chr21', 'chr22']):
        for i, chrom in enumerate(['chr13', 'chr14', 'chr15', 'chr21', 'chr22']):
            # total_snvs[chrom] = {'father': [0, 0, 0, 0, 0], 'mother': [0, 0, 0, 0, 0]}
            for g2 in ['mother', 'father']:
                # snv_path = f'{workdir}/snv_calls/{sample}_from_father.{g2}.snvs.tsv'
                valid_path = f'{var_dir}/{aligner}_val/validated_calls/{sample}_from_{g2}.{chrom}.validated_snvs.tsv'

                if not os.path.exists(valid_path):
                    continue
                # raw_snv_tbl = pd.read_csv(snv_path, sep='\t')
                valid_tbl = pd.read_csv(valid_path, sep='\t', header=[0])
                valid_tbl = valid_tbl.loc[valid_tbl['final_validation']=='true_de_novo']

                if len(valid_tbl) == 0:
                    continue

                ref_tig = valid_tbl['chr'].values[0]
                # ref_qarm_start = qarm_range.at[ref_tig, 'start']
                ref_qarm_start = int(qarm_range.at[ref_tig, 'qARM_START'])
                ref_qarm_side_base = qarm_range.at[ref_tig, 'size']
                ref_qarm_asat_base = qarm_range.at[ref_tig, 'asat_size']

                ## Exclude DNMs on contigs with limited p-arm sequences
                if ref_qarm_side_base / contig_size.at[ref_tig, 'size'] > 0.8:
                    continue

                qry_tig = g3_tig_matches.loc[g3_tig_matches['G2-TIG'].str.contains(ref_tig)]['G3-TIG'].values[0]
                vcf_tbl = pd.read_csv(f'{var_dir}/{qry_tig}/{qry_tig}.{aligner}.tsv', sep='\t', index_col=['id'])


                aln_tbl = pd.read_csv(f'./classify_alns/qry_contig_pafs/{qry_tig}.iden99.aln1Mb.tsv', sep='\t', header=0)
                aln_tbl['block_size'] = aln_tbl.apply(lambda row: row['ref_end'] - row['ref_start'], axis=1)

                slider_failed = pd.read_csv(f'{var_dir}/{qry_tig}/slider_align_w10k.failed.wfmash.bed', sep='\t', names=['ref', 'start', 'end'])
                failed_ref_dict = {}
                if len(slider_failed) > 0:
                    for idx, row in slider_failed.iterrows():
                        if row['ref'] not in failed_ref_dict:
                            failed_ref_dict[row['ref']] = IntervalTree()
                        failed_ref_dict[row['ref']][row['start']: row['end']] = [row['start'], row['end']]

                valid_tbl = valid_tbl.loc[valid_tbl['final_validation'] == 'true_de_novo']
                ref_chrom = ref_tig.split('_')[-1]

                for _, row in aln_tbl.iterrows():
                    callable_stats.append([row['ref'], row['ref_start'], row['ref_end'], qry_tig, 'transmission'])
                    bad_regions = list(flagged_dict[ref_tig].overlap(row['ref_start'], row['ref_end'])) if ref_tig in flagged_dict else []
                    bad_bases = 0
                    if bad_regions:
                        for itvl in bad_regions:
                            bad_bases += get_overlaps(row['ref_start'], row['ref_end'], itvl.begin, itvl.end)
                            callable_stats.append([row['ref'], itvl.begin, itvl.end, qry_tig, 'masked'])

                    if failed_ref_dict and ref_tig in failed_ref_dict:
                        slider_regions = list(failed_ref_dict[ref_tig].overlap(row['ref_start'], row['ref_end']))
                        if slider_regions:
                            for itvl in slider_regions:
                                bad_bases += get_overlaps(row['ref_start'], row['ref_end'], itvl.begin, itvl.end)
                                callable_stats.append([row['ref'], itvl.begin, itvl.end, qry_tig, 'slider_failed'])

                    bases[g2] += row['block_size'] - bad_bases
                    bases_ex_aSat[g2] += row['block_size'] - bad_bases
                    bases_chrom[g2][ref_chrom] += row['block_size'] - bad_bases

                bases[g2] -= ref_qarm_side_base
                bases_ex_aSat[g2] -= ref_qarm_asat_base
                bases_chrom[g2][ref_chrom] -= ref_qarm_side_base

                callable_stats.append([ref_tig, ref_qarm_start, qarm_range.at[ref_tig, 'end'], qry_tig, 'qarm'])

                pass_region_dnms = []

                for idx, row in valid_tbl.iterrows():
                    ref_pos = row['pos']
                    ## variation on qarm side
                    if ref_pos > ref_qarm_start:
                        continue
                    ## Variation in nucfreq or flagger masked regions
                    if ref_tig in flagged_dict and flagged_dict[ref_tig].overlaps(ref_pos, ref_pos+1):
                        continue
                    ## Variation in slider align failed regions
                    if ref_tig in failed_ref_dict and failed_ref_dict[ref_tig].overlaps(ref_pos, ref_pos + 1):
                        continue
                    qry_pos = vcf_tbl.at[row['id'], 'qry_pos']
                    ## Variation in nucfreq or flagger masked regions
                    if qry_tig in flagged_dict and flagged_dict[qry_tig].overlaps(qry_pos, qry_pos + 1):
                        continue
                    ref_hap = '_'.join(row['chr'].split('_')[0:2])

                    de_novo[g2][i] += 1
                    pass_region_dnms.append(row.tolist())
                    if g2 == 'mother':
                        mother_dnms.append([row['chr'], row['pos'], row['pos'] + 1, row['ref'], row['alt'], row['id'], qry_tig, sample, ref_hap, row['chr'].split('_')[-1]])
                    else:
                        fater_dnms.append([row['chr'], row['pos'], row['pos'] + 1, row['ref'], row['alt'], row['id'], qry_tig, sample, ref_hap, row['chr'].split('_')[-1]])

        for i, chrom in enumerate(['chr13', 'chr14', 'chr15', 'chr21', 'chr22']):
            chrom_pat_bases = bases_chrom['father'][chrom]
            chrom_mat_bases = bases_chrom['mother'][chrom]
            pat_rate = 0 if chrom_pat_bases == 0 else de_novo['father'][i]/chrom_pat_bases
            mat_rate = 0 if chrom_mat_bases == 0 else de_novo['mother'][i]/chrom_mat_bases
            mut_chrom.append([de_novo['father'][i], chrom_pat_bases, pat_rate, chrom, sample, 'Pat'])
            mut_chrom.append([de_novo['mother'][i], chrom_mat_bases, mat_rate, chrom, sample, 'Mat'])

        pat_bases = bases['father']
        mat_bases = bases['mother']

        print(f'-----{sample}------')
        print('\tPaternal SNV:', sum(de_novo['father']))
        print('\tMaternal SNV:', sum(de_novo['mother']))
        print('\tPaternal bases:', pat_bases)
        print('\tMaternal bases:', mat_bases)

        print('\tPaternal rate:', sum(de_novo['father']) / bases['father'])
        print('\tMaternal rate:', sum(de_novo['mother']) / bases['mother'])

        mut_rate.append([sum(de_novo['father']), bases['father'], bases_ex_aSat['father'], sum(de_novo['father']) / pat_bases, sample, 'Pat'])
        mut_rate.append([sum(de_novo['mother']), bases['mother'], bases_ex_aSat['mother'], sum(de_novo['mother']) / mat_bases, sample, 'Mat'])

    df_rate = pd.DataFrame(mut_rate, columns=['count', 'bases', 'bases_ex_asat', 'rate', 'g3_sample', 'hap'])
    df_rate.to_csv(f'{var_dir}/{aligner}_val/validated_calls/dnms_rate.tsv', sep='\t', header=True,index=False)

    df_chrom = pd.DataFrame(mut_chrom, columns=['count', 'bases', 'rate', 'chrom', 'g3_sample', 'hap'])
    df_chrom.to_csv(f'{var_dir}/{aligner}_val/validated_calls/dnms_rate_bychrom.tsv', sep='\t', header=True,index=False)

    all_dnm_df = pd.DataFrame(mother_dnms + fater_dnms, columns=['ref', 'pos', 'end', 'ref', 'alt', 'id', 'qry_tig', 'qry_sample', 'ref_hap','chrom'])
    all_dnm_df.to_csv(f'{var_dir}/{aligner}_val/validated_calls/all_dnms.tsv', sep='\t', header=True,index=False)

def get_valid_sv(var_dir):
    qarm_range = pd.read_csv('/annotation/contig_qc_summary.tsv', sep='\t', names=['tig_name', 'start', 'end'])
    qarm_range['size'] = qarm_range.apply(lambda row:row['end'] - row['start'], axis=1)
    qarm_range.set_index('tig_name', inplace=True)

    flagged_dict = {}
    for line in open('./annotation/flagger_nucfreq_merged.bed'):
        chrom, start, end = line.strip().split('\t')
        if chrom not in flagged_dict:
            flagged_dict[chrom] = IntervalTree()
        flagged_dict[chrom][int(start): int(end)] = [int(start), int(end)]

    mother_dn, father_dn = [], []

    validate_tbl = pd.read_csv(f'{var_dir}/validate_sv.tab', sep='\t', header=0)
    for idx, row in validate_tbl.iterrows():
        qry_tig = row['qry_tig']
        sample, chrom, parent = row['sample'], row['chrom'], row['parent']
        valid_res = pd.read_csv(f'{var_dir}/results/{sample}_{chrom}_{parent}/validate_res.tsv', sep='\t', header=0)
        denovos = valid_res.loc[valid_res['final_validation'] == 'true_de_novo']
        if len(denovos) == 0:
            continue

        slider_failed = pd.read_csv(f'variants/{qry_tig}/slider_align_w10k.failed.wfmash.bed', sep='\t', names=['ref', 'start', 'end'])
        failed_ref_dict = {}
        if len(slider_failed) > 0:
            for idx, row in slider_failed.iterrows():
                if row['ref'] not in failed_ref_dict:
                    failed_ref_dict[row['ref']] = IntervalTree()
                failed_ref_dict[row['ref']][row['start']: row['end']] = [row['start'], row['end']]

        for _, drow in denovos.iterrows():
            ref_tig = drow['chrom']

            ref_qarm_start = qarm_range.at[ref_tig, 'start']
            ref_start, ref_end = drow['start'], drow['end']

            ## variation on qarm side
            if ref_start > ref_qarm_start:
                continue

            ## Variation in nucfreq or flagger masked regions
            if ref_tig in flagged_dict and flagged_dict[ref_tig].overlaps(ref_start, ref_end):
                continue
            ## Variation in slider align failed regions
            if ref_tig in failed_ref_dict and failed_ref_dict[ref_tig].overlaps(ref_start, ref_end):
                continue

            ## Variation in flagged regions on query contig
            if qry_tig in flagged_dict and flagged_dict[qry_tig].overlaps(ref_start, ref_end):
                continue

            if parent == 'mother':
                mother_dn.append(drow.tolist() + [ref_tig, qry_tig])
            else:
                father_dn.append(drow.tolist() + [ref_tig, qry_tig])

    all_dnm_df = pd.DataFrame(mother_dn + father_dn,columns=['chrom', 'start', 'end', 'length', 'svtype', 'final_validation', 'supp_info', 'ref_tig', 'qry_tig'])
    all_dnm_df.to_csv(f'{var_dir}/all_dnsvs.tsv', sep='\t', header=True, index=False)

def main():
    g3_samples = ['NA12879', 'NA12881', 'NA12882', 'NA12883', 'NA12884','NA12885', 'NA12886', 'NA12887']
    g4_samples = ['K200081', 'K200082', 'K200084', 'K200085', 'K200087']

    ## Incorrect alignment
    g2g3_false_trans = ['NA12887_1_haplotype1-0000013_chr13', 'NA12882_2_haplotype2-0000107_chr21', 'NA12882_2_haplotype2-0000112_chr14']
    g3_NA12879_false_recomb = ['K200086_1_haplotype1-0000033_chr22', 'K200082_1_haplotype1-0000020_chr21']

    ## Get HiFi and ONT reads validated de novo SNV
    get_valid_dnms('variants_dir/', 'wfmash', g3_samples, g2g3_false_trans)
    get_valid_dnms('variants_dir/', 'wfmash', g4_samples, g3_NA12879_false_recomb)

    ## Get HiFi and ONT reads validated de novo SV
    get_valid_sv(f'variants_dir/sv_val')

if __name__ == '__main__':
    main()
