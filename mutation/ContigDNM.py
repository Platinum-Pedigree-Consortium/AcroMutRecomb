#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/4/25
'''
import pandas as pd
from intervaltree import IntervalTree


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


if __name__ == '__main__':
    main()
