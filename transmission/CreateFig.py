#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/6/25
'''

import json
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import seaborn as sns
from matplotlib.lines import Line2D

G2 = ['NA12877', 'NA12878']
G3 = ['NA12879', 'NA12881', 'NA12882', 'NA12883', 'NA12884','NA12885', 'NA12886', 'NA12887', 'K200080', 'K200100']
G4 = ['K200081', 'K200082', 'K200084', 'K200085', 'K200086', 'K200087', 'K200101', 'K200102', 'K200103', 'K200104', 'K200106']

SAT_COLORS = {'bSat': '#FFD63A', 'HSat1A': '#00DE60', 'HSat1B': '#1B998B','HSat3': '#fdb462', 'ACRO': '#3E7C17', 'aSat': '#990000', 'PHR_LOC': '#FF6600',
              'Unlabelled': '#bdbdbd', 'Others': '#BED7DC', 'DJ': '#990099', 'DJflank': '#990099', 'DJarm': '#990099', 'PJ': '#FF6666',
              'SegDup': '#386cb0', 'rDNA': '#9966FF', 'SST1': '#FF2DD1', 'parm-sd': '#386cb0', 'qarm-sd': '#386cb0', 'SEQ': '#f0f0f0',
              'SD': '#386cb0', 'aSat_HOR': '#9970ab', 'CER': '#00CCCC', 'GAP': '#252525'}


def figure3a():
    contig_size = pd.read_csv('./all_pq_contigs.size', sep='\t', usecols=[0,1], names=['contig', 'size'])


    g3_child = ['NA12879_1', 'NA12879_2', 'NA12881_1', 'NA12881_2',
                'NA12882_1', 'NA12882_2', 'NA12883_1', 'NA12883_2',
                'NA12884_1', 'NA12884_2', 'NA12885_1', 'NA12885_2',
                'NA12886_1', 'NA12886_2', 'NA12887_1', 'NA12887_2']
    NA12879_g4_child = ['K200081_1', 'K200081_2', 'K200082_1', 'K200082_2', 'K200084_1', 'K200084_2', 'K200085_1',
                        'K200085_2', 'K200087_1', 'K200087_2']

    child_haps = g3_child + NA12879_g4_child
    child_name = []
    for ele in child_haps:
        sample = ele.split('_')[0]
        if sample not in child_name:
            child_name.append(sample)

    samples_bases = {'Total': [0 for _ in range(len(child_haps))], '1to1': [0 for _ in range(len(child_haps))],
                     '1toN': [0 for _ in range(len(child_haps))], 'Others': [0 for _ in range(len(child_haps))]}

    sample_parm_bases = [0 for _ in range(len(child_haps))]

    for idx, row in contig_size.iterrows():
        sample_hap = row['contig'].split('_')[0] + '_' + row['contig'].split('_')[1]
        if sample_hap not in child_haps:
            continue
        samples_bases['Total'][child_haps.index(sample_hap)] += row['size']/1000000

    g3_trans_dict = json.load(open(f'./wfmash_g2g3/aligned_bases.json', 'r'))
    g3_parm_trans_dict = json.load(open(f'./wfmash_g2g3/aligned_bases_parm.json', 'r'))
    NA12879_g4_trans_dict = json.load(open(f'./wfmash_g3_NA12879/aligned_bases.json', 'r'))
    NA12879_g4_parm_trans_dict = json.load(open(f'./wfmash_g3_NA12879/aligned_bases_parm.json', 'r'))

    ## get all parm, qarm transmitted bases
    last_idx = 0
    for aln_type, bases_list in g3_trans_dict.items():
        if aln_type == 'Total':
            continue
        for i, base in enumerate(bases_list):
            samples_bases[aln_type][i] += base
            last_idx = i

    last_idx += 1
    for aln_type, bases_list in NA12879_g4_trans_dict.items():
        new_base_list = []
        for ele in bases_list:
            if ele > 30:
                new_base_list.append(ele)

        if aln_type == 'Total':
            continue
        for j, base in enumerate(new_base_list):
            samples_bases[aln_type][last_idx + j] += base

    ## get the parm transmitted bases
    last_idx = 0
    for aln_type, bases_list in g3_parm_trans_dict.items():
        if aln_type == 'Total':
            continue
        for i, base in enumerate(bases_list):
            sample_parm_bases[i] += base
            last_idx = i

    last_idx += 1
    for aln_type, bases_list in NA12879_g4_parm_trans_dict.items():
        new_base_list = []
        for ele in bases_list:
            if ele > 0:
                new_base_list.append(ele)
        if aln_type == 'Total':
            continue
        for j, base in enumerate(new_base_list):
            sample_parm_bases[last_idx + j] += base


    g3_bases = sum(samples_bases['1to1'][0:len(g3_child)])
    g4_NA12879_bases = sum(samples_bases['1to1'][len(g3_child)+1:])

    print(g3_bases)
    print(g4_NA12879_bases)

    print(sum(sample_parm_bases[0: len(g3_child)]))
    print(sum(sample_parm_bases[len(g3_child)+1:]))
    print(g4_NA12879_bases)
    print("Median per-hap:", np.median(samples_bases['1to1']))
    print("Median parm per-hap:", np.median(sample_parm_bases))

    legends = [Patch(label='Transmitted', color='#fd8d3c'), Patch(label='Others', color='#74c476')]

    xticks = np.arange(len(child_haps))
    width = 0.6
    fig, ax = plt.subplots(figsize=(8, 7))
    # other_bases = [i - j - k for i, j, k in zip(samples_bases['Total'], samples_bases['1to1'], samples_bases['1toN'])]
    ax.barh(xticks, samples_bases['Total'], height=width, color='#74c476', label='Total')
    # ax.bar(xticks, samples_bases['Others'], width=width, color='#6baed6', label='Others')
    ax.barh(xticks, samples_bases['1to1'], height=width, color='#fd8d3c', label='1to1')
    ax.barh(xticks, sample_parm_bases, height=width, color='#4393c3', label='parm')
    # ax.bar(xticks, samples_bases['1toN'],
    #        bottom=[i + j for i, j in zip(samples_bases['1to1'], samples_bases['Others'])], width=width, color='#32974c',
    #        label='1toN')

    # ax.legend(loc='lower left')
    ax.legend(handles=legends, ncols=2)
    ax.set_xticks([0, 20, 40, 60, 80], labels=[0, 20, 40, 60, 80], fontsize=12)
    # ax.set_yticklabels([0, 2000, 4000, 6000])
    ax.set_xlabel('Transmitted bases (Mbp)', fontsize=13)
    yticks = []
    sample_ticks = np.arange(len(child_name))
    for i, ele in enumerate(child_name):
        yticks.append((sample_ticks[i] + (i + 1)*width))
    ax.set_yticks(yticks)
    ax.set_yticklabels(child_name, fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    plt.show()

def figure4a():
    chrom_val = {'chr13': 1, 'chr14': 2, 'chr15': 3, 'chr21': 4, 'chr22': 5}
    rec_tbl = pd.read_csv('./recomb_info.tsv', sep='\t', header=0)
    qc_tbl = pd.read_csv('./contig_qc_summary.tsv', sep='\t', index_col=['HAP'])
    contig_size = pd.read_csv('./all_pq_contigs.size', sep='\t', names=['contig', 'size'], usecols=[0, 1])
    contig_size.set_index('contig', inplace=True)

    qarm_rec = rec_tbl.loc[(rec_tbl['Arm']=='qarm')|(rec_tbl['Arm']=='parm')]
    qarm_rec['dist2aSat'] = qarm_rec.apply(lambda row: (row['child_contig_pos'] - int(qc_tbl.at[row['Child_contig'], 'qARM_START']))/1000000, axis=1)
    qarm_rec['Tel2aSat'] = qarm_rec.apply(lambda row: (1 - int(qc_tbl.at[row['Child_contig'], 'qARM_START']))/1000000, axis=1)
    qarm_rec['End2aSat'] = qarm_rec.apply(lambda row: (contig_size.at[row['Child_contig'], 'size'] - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) / 1000000, axis=1)
    qarm_rec['distal'] = qarm_rec.apply(lambda row: qc_tbl.at[row['Child_contig'], 'DISTAL'], axis=1)
    qarm_rec['val'] = qarm_rec.apply(lambda row: chrom_val[row['Chrom']], axis=1)


    qarm_has_distal = qarm_rec.loc[qarm_rec['distal']=='Yes']
    qarm_has_distal['dj2aSat'] = qarm_has_distal.apply(lambda row: (int(qc_tbl.at[row['Child_contig'], 'DJ_END']) - int(qc_tbl.at[row['Child_contig'], 'qARM_START']))/1000000, axis=1)


    prdm9 = pd.read_csv('./all_recomb_seq/fimo_out/best_site.narrowPeak', sep='\t', usecols=[0,1,2], names=['contig', 'start', 'end'])

    fig, (ax_p, ax) = plt.subplots(2, 1, sharex=True, figsize=(13, 6), gridspec_kw={'height_ratios': [1, 7]})
    # parm_rc_contigs = ['NA12879_1_haplotype1-0000002_chr13', 'K200084_1_haplotype1-0000012_chr21', 'NA12879_1_haplotype1-0000019_chr21']
    parm_rc_info = rec_tbl.loc[rec_tbl['Child_contig']=='K200084_1_haplotype1-0000012_chr21']


    g4_tel2asat = (1 - int(qc_tbl.at['K200084_1_haplotype1-0000012_chr21', 'qARM_START']))/1000000
    g4_dist2asat = (parm_rc_info['child_contig_pos'].values[0] - int(qc_tbl.at['K200084_1_haplotype1-0000012_chr21', 'qARM_START']))/1000000
    g4_End2aSat = (contig_size.at['K200084_1_haplotype1-0000012_chr21', 'size'] - int(qc_tbl.at['K200084_1_haplotype1-0000012_chr21', 'qARM_START'])) / 1000000

    g3_chr13_tel2asat = (1 - int(qc_tbl.at['NA12879_1_haplotype1-0000002_chr13', 'qARM_START'])) / 1000000
    g3_chr13_dist2asat = (parm_rc_info['parent_contig_1_pos'].values[0] - int(qc_tbl.at['NA12879_1_haplotype1-0000002_chr13', 'qARM_START'])) / 1000000
    g3_chr13_End2aSat = (contig_size.at['NA12879_1_haplotype1-0000002_chr13', 'size'] - int(qc_tbl.at['NA12879_1_haplotype1-0000002_chr13', 'qARM_START'])) / 1000000

    g3_chr21_tel2asat = (1 - int(qc_tbl.at['NA12879_1_haplotype1-0000019_chr21', 'qARM_START'])) / 1000000
    g3_chr21_dist2asat = (parm_rc_info['parent_contig_2_pos'].values[0] - int(qc_tbl.at['NA12879_1_haplotype1-0000019_chr21', 'qARM_START'])) / 1000000
    g3_chr21_End2aSat = (contig_size.at['NA12879_1_haplotype1-0000019_chr21', 'size'] - int(qc_tbl.at['NA12879_1_haplotype1-0000019_chr21', 'qARM_START'])) / 1000000

    contig_peak_tbl = prdm9.loc[prdm9['contig'] == 'K200084_1_haplotype1-0000012_chr21']
    peaks = []
    for _, peak_row in contig_peak_tbl.iterrows():
        peak_pos = (peak_row['start'] - int(qc_tbl.at['K200084_1_haplotype1-0000012_chr21', 'qARM_START'])) / 1000000
        dist2bp = g4_dist2asat - peak_pos
        peaks.append([abs(dist2bp), peak_pos])


    ax_p.plot([g3_chr13_tel2asat, g3_chr13_End2aSat], [0.3, 0.3], color='#4d4d4d', lw=1)
    ax_p.scatter([g3_chr13_tel2asat], [0.3], marker='o', color='b')
    ax_p.scatter([g3_chr13_dist2asat], [0.3], marker='^', color='r')
    ax_p.scatter([g3_chr13_End2aSat], [0.3], marker='o', color='b')

    ax_p.plot([g4_tel2asat, g4_End2aSat], [0.2, 0.2], color='#4d4d4d', lw=1)
    ax_p.scatter([g4_tel2asat], [0.2], marker='o', color='b')
    ax_p.scatter([g4_dist2asat], [0.2], marker='^', color='r')
    ax_p.scatter([g4_End2aSat], [0.2], marker='o', color='b')

    sorted_peaks = sorted(peaks, key=lambda x: x[0])
    ax_p.scatter(sorted_peaks[0][1], [0.2], marker='|', color='k', s=50)

    ax_p.plot([g3_chr21_tel2asat, g3_chr21_End2aSat], [0.1, 0.1], color='#4d4d4d', lw=1)
    ax_p.scatter([g3_chr21_tel2asat], [0.1], marker='o', color='b')
    ax_p.scatter([g3_chr21_dist2asat], [0.1], marker='^', color='r')
    ax_p.scatter([g3_chr21_End2aSat], [0.1], marker='o', color='b')
    ax_p.spines[['left', 'top', 'right', 'bottom']].set_visible(False)
    ax_p.set_yticks([0.1,0.2,0.3], labels=['NA12879_1', '200084_1', 'NA12879_1'], fontsize=13)
    ax_p.axvline(x=0, ls='--', lw=1.5, color='grey')
    ax_p.set_ylim([0.05, 0.35])
    all_peaks = []
    ylabels = []
    yticks = []
    yvals = 1
    qarm_rec.sort_values('Chrom', inplace=True)

    for idx, row in qarm_rec.iterrows():
        if row['Child_contig'] == 'K200084_1_haplotype1-0000012_chr21':
            continue
        child_hap = row['Child_contig'].split('_')[0] + '_' + row['Child_contig'].split('_')[1]
        label = child_hap + '_' + row['Chrom']
        ylabels.append(label)
        ax.plot([row['Tel2aSat'], row['End2aSat']], [yvals, yvals], color='#4d4d4d', lw=1)

        ax.scatter([row['Tel2aSat']], [yvals], marker='o', color='b')
        ax.scatter([row['dist2aSat']], [yvals], marker='^', color='r')
        ax.scatter([row['End2aSat']], [yvals], marker='o', color='b')

        contig_peak_tbl = prdm9.loc[prdm9['contig']==row['Child_contig']]
        peaks = []
        for _, peak_row in contig_peak_tbl.iterrows():
            peak_pos = (peak_row['start'] - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) / 1000000
            dist2bp = row['dist2aSat'] - peak_pos
            peaks.append([abs(dist2bp), peak_pos, abs(row['child_contig_pos'] - peak_row['start'])/1000])
            all_peaks.append([abs(row['child_contig_pos'] - peak_row['start'])/1000000, peak_pos, row['Child_contig'], row['Child_contig'].split('_')[-1]])

        sorted_peaks = sorted(peaks, key=lambda x:x[0])
        print(row['Child_contig'], sorted_peaks[0][2])
        ax.scatter(sorted_peaks[0][1], [yvals], marker='|', color='k', s=50)

        ## Add inversion breakpoint to three haplotypes
        if row['Child_contig'] == 'NA12886_1_haplotype1-0000027_chr15':
            ax.scatter((11192184 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) / 1000000, [yvals], marker='d', color='g')
            ax.scatter((12498225 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) /1000000, [yvals], marker='d', color='g')

        if row['Child_contig'] == 'NA12887_2_haplotype2-0000070_chr15':
            ax.scatter((13083530 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) / 1000000, [yvals], marker='d', color='g')
            ax.scatter((14401853 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) /1000000, [yvals], marker='d', color='g')

        if row['Child_contig'] == 'K200087_2_haplotype2-0000206_chr15':
            ax.scatter((8900438 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) / 1000000, [yvals], marker='d', color='g')
            ax.scatter((10194279 - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])) /1000000, [yvals], marker='d', color='g')

        # if row['Child_contig'] in inv_qry2:
        #     qry_inv_start = inv_tbl.loc[inv_tbl['qry2'] == row['Child_contig']]['qry2_start'].values[0] - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])
        #     qry_inv_end = inv_tbl.loc[inv_tbl['qry2'] == row['Child_contig']]['qry2_end'].values[0] - int(qc_tbl.at[row['Child_contig'], 'qARM_START'])
        #     ax.scatter(qry_inv_start/1000000, [yvals], marker='d', color='g')
        #     ax.scatter(qry_inv_end/1000000, [yvals], marker='d', color='g')

        yticks.append(yvals)
        yvals += 0.1

    ax.axvline(x=0, ls='--', lw=1.5, color='grey')
    # ax.set_xticks([-12, -8, -4, 0, 2, 4, 6], labels=[-10, -8, -6, -4, -2, 0, 2, 4, 6], fontsize=13)
    ax.set_ylabel('')
    ax.set_yticks(yticks, labels=[ele.split('_')[0] + '_' + ele.split('_')[1] for ele in ylabels], fontsize=13)
    ax.set_xlabel('Distance to aSat end (Mbp)', fontsize=14)
    ax.spines[['left', 'top', 'right']].set_visible(False)
    ax.legend('', frameon=False)
    fig.tight_layout()
    plt.show()

def figure5b():
    flagger_tbl = pd.read_csv('./contig_qc/pq_flagger_regions.bed', sep='\t',
                              usecols=[0, 1, 2], names=['contig', 'start', 'stop'])

    nucfreq_tbl = pd.read_csv('./contig_qc/pq_nucfreq_regions.bed', sep='\t',
                              usecols=[0, 1, 2], names=['contig', 'start', 'stop'])

    rec_tbl = pd.read_csv('./recomb/recomb_info.tsv', sep='\t', header=0)
    contig_size = pd.read_csv('./all_pq_contigs.size', sep='\t',
                              usecols=[0, 1], names=['contig', 'size'])
    contig_size.set_index('contig', inplace=True)

    contig = 'K200084_1_haplotype1-0000012_chr21'
    test_tble = rec_tbl.loc[rec_tbl['Child_contig']==contig]
    censat_tbl = pd.read_csv('./all_samples_censat.bed', sep='\t',
                             names=['contig', 'start', 'end', 'label', 'color'])

    for idx, row in test_tble.iterrows():
        child_contig = row['Child_contig']
        sample = child_contig.split('_')[0]

        wr_path = './contig_aln/wfmash_g2g3'
        if sample in G4:
            wr_path = './contig_aln/wfmash_g3_NA12879'

        qry_ref_matches = pd.read_csv(f'{wr_path}/validate_aln.tab', sep='\t', index_col=['QUERY'])

        ref_contig = qry_ref_matches.at[child_contig, 'TARGET']

        tmp_contig = ref_contig.split('_')

        ref_contig_1, ref_contig_2 = f'{tmp_contig[0]}_{tmp_contig[1]}_{tmp_contig[2]}_{tmp_contig[3]}', f'{tmp_contig[4]}_{tmp_contig[5]}_{tmp_contig[6]}_{tmp_contig[7]}'

        ref_contig_pos1 = row['parent_contig_1_pos'] if ref_contig_1 == row['parent_contig_1'] else row['parent_contig_2_pos']
        ref_contig_pos2 = row['parent_contig_2_pos'] if ref_contig_2 == row['parent_contig_2'] else row['parent_contig_1_pos']

        hifi_bg = pd.read_csv(f'{wr_path}/variants/{child_contig}/{ref_contig}.HiFi.bg', sep='\t', names=['contig', 'start', 'end', 'cov'])
        hifi_bg['start_mb'] = hifi_bg.apply(lambda row: row['start']/1000000, axis=1)
        hifi_bg['cov_new'] = hifi_bg.apply(lambda row: row['cov'] + math.log(800), axis=1)

        ont_bg = pd.read_csv(f'{wr_path}/variants/{child_contig}/{ref_contig}.ONT.bg', sep='\t',names=['contig', 'start', 'end', 'cov'])
        ont_bg['start_mb'] = ont_bg.apply(lambda row: row['start'] / 1000000, axis=1)
        ont_bg['cov_new'] = ont_bg.apply(lambda row: row['cov'] + math.log(800), axis=1)

        fig, axes = plt.subplots(2, 1, sharey=True, figsize=(14, 5))
        fig.suptitle(f'Child contig: {child_contig}', fontsize=14, y=0.97)

        for i, ele in enumerate([ref_contig_1, ref_contig_2]):
            size = contig_size.at[ele, 'size']

            this_fagger = flagger_tbl.loc[flagger_tbl['contig'] == contig]
            this_nucfreq = nucfreq_tbl.loc[nucfreq_tbl['contig'] == contig]
            flagger_ranges = []
            nucfreq_ranges = []
            bad_kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}
            for idx, row in this_fagger.iterrows():
                this_start = row['start'] / 1000000
                group = 'Failed'
                if row['stop'] <= size:
                    this_end = (row['stop'] - row['start']) / 1000000
                    this_color = '#f0f0f0' if group == 'passed' else "#878787"
                    flagger_ranges.append((this_start, this_end))
                    bad_kwargs['edgecolor'].append(this_color)
                    bad_kwargs['facecolor'].append(this_color)
            for idx, row in this_nucfreq.iterrows():
                this_start = row['start'] / 1000000
                group = 'Failed'
                if row['stop'] <= size:
                    this_end = (row['stop'] - row['start']) / 1000000
                    this_color = '#f0f0f0' if group == 'passed' else "#878787"
                    nucfreq_ranges.append((this_start, this_end))
                    bad_kwargs['edgecolor'].append(this_color)
                    bad_kwargs['facecolor'].append(this_color)

            axes[i].broken_barh([(0, size/1000000)], (0, math.log(400)), edgecolor='#f0f0f0', facecolor='#f0f0f0')
            axes[i].broken_barh(flagger_ranges, (0, math.log(400)), **bad_kwargs)

            axes[i].broken_barh([(0, size/1000000)], (math.log(400), math.log(200)), edgecolor='#f0f0f0', facecolor='#f0f0f0')
            axes[i].broken_barh(nucfreq_ranges, (math.log(400), math.log(200)), **bad_kwargs)

            this_tig_asat = censat_tbl.loc[censat_tbl['contig'] == ele]
            sat_ranges = []
            sat_kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}

            for idx, row in this_tig_asat.iterrows():
                this_start = row['start'] / 1000000
                if row['end'] <= size:
                    if row['label'] in ['gray1', 'green2', 'yellow1', 'red2']:
                        continue
                    this_end = (row['end'] - row['start']) / 1000000
                    this_color = row['color'] if row['label'] != 'SEQ' else '#f0f0f0'
                    if row['label'] in SAT_COLORS:
                        this_color = SAT_COLORS[row['label']]
                    if ['HSat1', 'HSat2', 'HSat3', 'bSat', 'aSat', 'SST1']:
                        sat_ranges.append((this_start, this_end))
                        sat_kwargs['edgecolor'].append(this_color)
                        sat_kwargs['facecolor'].append(this_color)

            axes[i].broken_barh([(0, size/1000000)], (0, math.log(600)), edgecolor='#f0f0f0', facecolor='#f0f0f0')
            axes[i].broken_barh(sat_ranges, (0, math.log(600)), **sat_kwargs)

        print(ref_contig_pos1, ref_contig_pos2)
        sns.lineplot(data=hifi_bg.loc[(hifi_bg['contig']==ref_contig_1)&(hifi_bg['start']<=ref_contig_pos1)], x='start_mb', y='cov_new', color='#a6cee3', ax=axes[0])
        sns.lineplot(data=hifi_bg.loc[(hifi_bg['contig'] == ref_contig_1) & (hifi_bg['start'] > ref_contig_pos1)], x='start_mb', y='cov_new', color='#1f78b4', ax=axes[0])

        sns.lineplot(data=ont_bg.loc[(ont_bg['contig'] == ref_contig_1) & (ont_bg['start'] <= ref_contig_pos1)],
                     x='start_mb', y='cov_new', color='#fdbf6f', ax=axes[0])
        sns.lineplot(data=ont_bg.loc[(ont_bg['contig'] == ref_contig_1) & (ont_bg['start'] > ref_contig_pos1)],
                     x='start_mb', y='cov_new', color='#ff7f00', ax=axes[0])

        axes[0].axvline(ref_contig_pos1/1000000, lw=2, ls='--', color='r')
        axes[0].set_xlabel(f'Position on {ref_contig_1} (Mbp)', fontsize=14)

        sns.lineplot(data=hifi_bg.loc[(hifi_bg['contig'] == ref_contig_2)&(hifi_bg['start']>=ref_contig_pos2)], x='start_mb', y='cov_new', color='#a6cee3', ax=axes[1])
        sns.lineplot(data=hifi_bg.loc[(hifi_bg['contig'] == ref_contig_2) & (hifi_bg['start'] <ref_contig_pos2)], x='start_mb', y='cov_new', color='#1f78b4', ax=axes[1])

        sns.lineplot(data=ont_bg.loc[(ont_bg['contig'] == ref_contig_2) & (ont_bg['start'] >= ref_contig_pos2)],
                     x='start_mb', y='cov_new', color='#fdbf6f', ax=axes[1])
        sns.lineplot(data=ont_bg.loc[(ont_bg['contig'] == ref_contig_2) & (ont_bg['start'] < ref_contig_pos2)],
                     x='start_mb', y='cov_new', color='#ff7f00', ax=axes[1])

        axes[1].axvline(ref_contig_pos2/1000000, lw=2, ls='--', color='r')
        axes[1].set_xlabel(f'Position on {ref_contig_2} (Mbp)', fontsize=14)

        legends = [Line2D([0], [0], ls='--', color='r', label='Breakpoint'), Line2D([0], [0], color='#1f78b4', label='Trans'), Line2D([0], [0], color='#a6cee3', label='NotTrans')]
        axes[0].legend(handles=legends, loc='upper left', ncols=3)

        for ax in axes:
            ax.margins(x=0.01, tight=True)
            ax.set_yscale('log')
            ax.set_ylabel('Coverage', fontsize=13)
            ax.spines[['right', 'top']].set_visible(False)
        fig.tight_layout()
    plt.show()

def figure5c():
    segdup = pd.read_csv('./contig_qc/segdup/all_contig_segdup.bed', sep='\t',
                         names=['contig', 'start', 'end'])
    rec_tbl = pd.read_csv('./contig_aln/recomb/recomb_info.tsv',
                          sep='\t', header=0)
    contig_size = pd.read_csv('./all_pq_contigs.size', sep='\t', usecols=[0, 1], names=['contig', 'size'])
    contig_size.set_index('contig', inplace=True)

    prdm_hit_density = pd.read_csv(
        './fimo_out/fimo_motif_1-14_2_20kb_bined_cnt.bed',
        sep='\t', names=['sequence_name', 'start', 'stop', 'value'])

    prdm_best = pd.read_csv(
        './fimo_out/best_site.narrowPeak',
        sep='\t', usecols=[0, 1, 2], names=['sequence_name', 'start', 'stop'])

    censat_tbl = pd.read_csv('./all_samples_censat.bed', sep='\t',
                             names=['contig', 'start', 'end', 'label', 'color'])

    # test_contig = rec_tbl.loc[rec_tbl['Child_contig']=='K200084_1_haplotype1-0000012_chr21']

    for idx, row in rec_tbl.iterrows():
        child_contig = row['Child_contig']
        size = contig_size.at[child_contig, 'size']
        this_contig_hits = prdm_hit_density.loc[prdm_hit_density['sequence_name'] == child_contig]

        this_contig_best = prdm_best.loc[prdm_best['sequence_name']==child_contig]
        this_contig_hits['position'] = this_contig_hits.apply(lambda row: row['start']/1000000, axis=1)
        max_value = this_contig_hits['value'].max()
        if max_value > 500:
            this_contig_hits['new_value'] = this_contig_hits.apply(lambda row: row['value']/10, axis=1)
        else:
            this_contig_hits['new_value'] = this_contig_hits.apply(lambda row: row['value'], axis=1)
        print(child_contig)
        child = child_contig.split('_')[0]
        aln_dir = './wfmash_g2g3/variants/{child_contig}/mm2_ava'
        if child in G4:
            aln_dir = './wfmash_g3_NA12879/variants/{child_contig}/mm2_ava'

        prdm9_xranges = []
        prdm9_kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}

        for idx, row in this_contig_best.iterrows():
            start = row['start']/1000000
            # if row['q-value']<1e-8:
            end = (row['stop'] - row['start'])/1000000
            prdm9_kwargs['edgecolor'].append('k')
            prdm9_kwargs['facecolor'].append('k')
            prdm9_xranges.append((start, end))

        this_contig_sd = segdup.loc[segdup['contig'] == child_contig]
        sd_ranges = []
        sd_kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}
        sd_start_tracker = set()
        sd_end_tracker = []
        for idx, row in this_contig_sd.iterrows():
            this_start = row['start']
            if this_start in sd_start_tracker:
                continue

            if row['end'] <= size:
                # sd_last = row['end']
                this_end = (row['end'] - this_start)/1000000
                this_color = SAT_COLORS['SegDup']
                sd_ranges.append((this_start/ 1000000, this_end))
                sd_kwargs['edgecolor'].append(this_color)
                sd_kwargs['facecolor'].append(this_color)
                sd_start_tracker.add(this_start)
                sd_end_tracker.append(row['end'])

        aln_tbl = pd.read_csv(f'{aln_dir}/child_to_parents_10K_slider.tsv', sep='\t',
                              names=['qry', 'qry_start', 'qry_end', 'ref', 'ref_start', 'ref_end', 'iden'])

        aln_tbl['length'] = aln_tbl.apply(lambda row: row['qry_end'] - row['qry_start'], axis=1)
        aln_tbl.sort_values('qry_start', inplace=True)
        aln_tbl = aln_tbl.loc[(aln_tbl['ref'] != child_contig) & (aln_tbl['qry'] == child_contig)]
        tmp_aln_tbl = aln_tbl.loc[aln_tbl['iden']==100]

        colors = ['#66a61e', '#d95f02']
        # colors = ['#8da0cb', '#fc8d62']
        xranges = {}
        kwargs = {}
        color_dict = {}
        for idx, row in aln_tbl.iterrows():
            ref_name = row['ref']
            if ref_name not in xranges:
                xranges[ref_name] = []

            if ref_name not in kwargs:
                kwargs[ref_name] = {'edgecolor': [], 'facecolor': []}
                color_dict[ref_name] = colors[len(kwargs) - 1]

            this_start = row['qry_start']/1000000
            this_end = (row['qry_end'] - row['qry_start'])/1000000
            if row['length'] == 10000 and row['iden'] == 100:
                xranges[ref_name].append((this_start, this_end))
                kwargs[ref_name]['edgecolor'].append(color_dict[ref_name])
                kwargs[ref_name]['facecolor'].append(color_dict[ref_name])

        qry_sat = censat_tbl.loc[censat_tbl['contig'] == child_contig]
        sat_ranges = []
        sat_kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}
        for idx, row in qry_sat.iterrows():
            this_start = row['start']/1000000
            if row['end'] <= size:
                if row['label'] not in ['HSat1A', 'HSat1B', 'HSat2', 'HSat3', 'bSat', 'aSat', 'SST1', 'rDNA', 'GAP']:
                    continue
                this_end = (row['end'] - row['start'])/1000000
                if row['label'] in SAT_COLORS:
                    this_color = SAT_COLORS[row['label']]
                    sat_ranges.append((this_start, this_end))
                    sat_kwargs['edgecolor'].append(this_color)
                    sat_kwargs['facecolor'].append(this_color)


        fig, ax = plt.subplots(figsize=(16, 4))
        ax.broken_barh([(0, size / 1000000)], (-30, -10), facecolor='#f0f0f0', edgecolor='#f0f0f0')
        ax.broken_barh(sd_ranges, (-30, -10), **sd_kwargs)

        ## Add PHR for chr13-chr21 ectopic recombination
        if child_contig == 'K200084_1_haplotype1-0000012_chr21':
            ax.broken_barh([(1505837/1000000, (4395789-1505837)/1000000)], (-30, -10), facecolor='#6baed6', edgecolor='#6baed6')

        ax.broken_barh([(0, size / 1000000)], (-20, -10), facecolor='#f0f0f0', edgecolor='#f0f0f0')
        ax.broken_barh(sat_ranges, (-20, -10), **sat_kwargs)

        ax.broken_barh([(0, size / 1000000)], (-10, -10), facecolor='#f0f0f0', edgecolor='#f0f0f0')
        for ref, aln_ranges in xranges.items():
            aln_kwargs = kwargs[ref]
            ax.broken_barh(aln_ranges, (-10, -10), **aln_kwargs)
        ax.broken_barh([(0, size / 1000000)], (0, -10), facecolor='#f0f0f0', edgecolor='#f0f0f0')
        ax.broken_barh(prdm9_xranges, (0, -10), **prdm9_kwargs)
        sns.lineplot(data=this_contig_hits, x='position', y='new_value', ax=ax)
        # ax.axvline(x=bp/1000000, color='r', lw=1)
        ax.spines[['right', 'top', 'left']].set_visible(False)
        ax.set_ylabel('PRDM9 motif hits', fontsize=14)
        ax.set_xlabel('Position (Mbp)', fontsize=14)
        ax.margins(x=0.01, y=0, tight=True)
        # ax.spines['bottom'].set_position(('data', 0))
        # ax.set_yticks([])
        fig.tight_layout()
    # plt.show()

def main():

    figure3a()
    figure4a()

    ## Create figure5b-c and all similar figures for other 18 recombination shown in Data S1
    figure5b()
    figure5c()

if __name__ == '__main__':
    main()
