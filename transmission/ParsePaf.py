#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/4/25
'''
import os
import re
import json
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
from intervaltree import IntervalTree


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

def insdel_from_cigar(cg_tuples):
    insdel = []
    snps = []
    readPos = 0
    refPos = 0
    for opLen, op in cg_tuples:
        if op == "N" or op == "S":
            readPos += opLen
        elif op == "I":
            insdel.append([refPos, readPos, opLen, 'INS'])
            readPos += opLen
        elif op == "D":
            insdel.append([refPos, readPos, opLen, 'DEL'])
            refPos += opLen
        elif op == 'X':
            # snps.append([refPos, readPos, 'SNP'])
            readPos += opLen
            refPos += opLen
        elif op in ["M", "E", '=']:
            readPos += opLen
            refPos += opLen
        elif op == 'H':
            pass
        else:
            pass
    return insdel

hap_order = ['NA12877_1', 'NA12877_2', 'NA12878_1', 'NA12878_2', 'NA12879_1', 'NA12879_2', 'NA12881_1',
             'NA12881_2','NA12882_1', 'NA12882_2', 'NA12883_1', 'NA12883_2', 'NA12884_1', 'NA12884_2', 'NA12885_1',
             'NA12885_2','NA12886_1', 'NA12886_2', 'NA12887_1', 'NA12887_2', 'K200080_1', 'K200080_2', 'K200081_1', 'K200081_2',
             'K200082_1', 'K200082_2', 'K200084_1', 'K200084_2', 'K200085_1', 'K200085_2', 'K200086_1', 'K200086_2', 'K200087_1', 'K200087_2',
             'K200100_1', 'K200100_2', 'K200101_1', 'K200101_2', 'K200102_1', 'K200102_2', 'K200103_1', 'K200103_2', 'K200104_1', 'K200104_2',
             'K200106_1', 'K200106_2']

def classify_aligns(acro_parms, child_haps, paf_dir, paf_path, aln_info, out_paf_dir, iden_thresh, min_aln_size, separate_single_qry=False):
    '''
    Extract alignments between G3 QRY and G2 REF:
    1. QRY aligned to a group of REFs, wfmash_wgs.qry_aligns.txt
    2. Best 1to1 alignment, one2one_dir/wfmash_wgs.paf
    '''
    qc_tbl = pd.read_csv('contig_qc_summary.tsv', sep='\t', index_col=['HAP'])

    ## Read all g2 and g3 acro-parm contigs ##
    acro_parm_tracker = {}
    hap_tigs = [0 for _ in range(len(hap_order))]
    for line in open(acro_parms):
        entries = line.strip().split('\t')
        acro_parm_tracker[entries[0]] = [int(entries[1]), 0]
        hap_name = entries[0].split('_')[0] + '_' + entries[0].split('_')[1]
        hap_tigs[hap_order.index(hap_name)] += int(entries[1]) / 1000000


    qry_tig_len, ref_tig_len = {}, {}
    qry_alns = {}
    qry_aln_blocks = []

    with open(paf_path, 'r') as fin:
        for paf_line in fin:
            # [rstart, rend, qstart, qend, rlen, qlen, iden, rdir, qdir, rchr, qchr, cg] = parse_paf_line(line)

            line = paf_line.strip().split()

            rchr = line[5]  ## G2 sequence name
            qchr = line[0]  ## G3 sequence name

            qry_tig_len[qchr] = int(line[1])
            ref_tig_len[rchr] = int(line[6])

            rstart = int(line[7]) + 1
            rend = int(line[8])

            qstart = int(line[2]) + 1
            qend = int(line[3])
            rlen = abs(rend - rstart) + 1
            qlen = abs(qend - qstart) + 1

            if qlen < min_aln_size:
                continue

            cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg'][0]

            cg_tuple = cigar_tuple(cg)
            iden = round((sum([int(i[0]) for i in cg_tuple if i[1] == '=' or i[1] == 'M']) / sum(
                [int(i[0]) for i in cg_tuple if i[1] in {'=', 'M', 'X', 'D', 'I'}])) * 100, 2)

            iden_tag = f'>={iden_thresh}' if iden >= iden_thresh else f'<{iden_thresh}'

            qry_aln_blocks.append([qlen, qchr.split('-')[1], iden_tag, iden])

            sco = round((qlen / int(line[1])) * iden, 2)
            if qchr not in qry_alns:
                qry_alns[qchr] = []

            qry_alns[qchr].append([[rstart, rend, qstart, qend, rlen, qlen, iden, sco, '+', line[4], rchr, qchr], paf_line.strip()])
            acro_parm_tracker[qchr][1] += 1
            acro_parm_tracker[rchr][1] += 1

    qry_contig_len = round(sum(qry_tig_len.values())/1000000, 2)
    ref_contig_len = round(sum(ref_tig_len.values())/1000000, 2)
    print(f'Total QRY contig length {qry_contig_len}Mb')
    print(f'Total REF contig length {ref_contig_len}Mb')


    ## Start to classify qry contigs into 1to1 and 1toN ##
    ######################################################

    aln_writer = open(aln_info, 'w')
    oneone_writer = open(f'{out_paf_dir}/iden{iden_thresh}.1to1.paf', 'w')
    onen_writer = open(f'{out_paf_dir}/iden{iden_thresh}.1ton.paf', 'w')
    # others_writer = open(f'{out_paf_dir}/g2g3_{region_tag}.iden{iden_thresh}.others.paf', 'w')

    all_qry_bases = 0
    oneone_bases = 0
    onen_bases = 0
    other_bases = 0
    bases_for_var = 0
    samples_contig = {'Total': [0 for _ in range(len(child_haps))], '1to1': [0 for _ in range(len(child_haps))], '1toN': [0 for _ in range(len(child_haps))]}
    samples_bases = {'Total': [0 for _ in range(len(child_haps))], '1to1': [0 for _ in range(len(child_haps))],
                     '1toN': [0 for _ in range(len(child_haps))], 'Others': [0 for _ in range(len(child_haps))]}

    samples_bases_parm = {'Total': [0 for _ in range(len(child_haps))], '1to1': [0 for _ in range(len(child_haps))],
                     '1toN': [0 for _ in range(len(child_haps))], 'Others': [0 for _ in range(len(child_haps))]}

    contig_sizes = []
    perfect_match_num = 0
    perfect_match_bases = 0
    perfect_contig = open(f'{out_paf_dir}/iden100.g3contig.txt', 'w')
    for qchr, alns in qry_alns.items():
        iden99_alns = []
        iden99_qry_bases = {}
        for ele in alns:
            if ele[0][6] >= iden_thresh:
                iden99_alns.append(ele)
                if ele[0][10] not in iden99_qry_bases:
                    iden99_qry_bases[ele[0][10]] = 0
                iden99_qry_bases[ele[0][10]] += ele[0][5]

        if len(iden99_alns) == 0:
            continue
        max_iden99_base = sorted(iden99_qry_bases.items(), key=lambda x:x[1], reverse=True)[0]

        ## Sort qry alignment by position
        sorted_alns = sorted(iden99_alns, key=lambda x:(x[0][2]), reverse=True)
        sorted_iden = sorted(iden99_alns, key=lambda x:(x[0][6]), reverse=True)
        if sorted_iden[0][0][6] == 100:
            qstart, qend = sorted_iden[0][0][2], sorted_iden[0][0][3]
            print(f'{qchr}\t{qstart}\t{qend}', file=perfect_contig)
            perfect_match_num += 1
            perfect_match_bases += sorted_iden[0][0][5]

        sorted_out = ''
        this_contig_aln_bases = 0
        sample_hap = qchr.split('_')[0] + '_' + qchr.split('_')[1]
        samples_contig['Total'][child_haps.index(sample_hap)] += 1

        if len(iden99_alns) == 1:
            if max_iden99_base[1] / qry_tig_len[qchr] < 0.5:
                continue
            this_aln, this_paf = sorted_alns[0]
            sorted_out = f'IDEN:{this_aln[6]},REF:{this_aln[10]},RSTART:{this_aln[0]},REND:{this_aln[1]},QSTART:{this_aln[2]},QEND:{this_aln[3]},ORI:{this_aln[8]}'
            print(f'{qchr}\t{acro_parm_tracker[qchr][0]}\t{len(sorted_alns)}\t{sorted_out}', file=aln_writer)
            if this_aln[6] >= iden_thresh:
                oneone_bases += this_aln[5]
                samples_bases['1to1'][child_haps.index(sample_hap)] += this_aln[5] / 1000000

                if qchr in qc_tbl.index:
                    samples_bases_parm['1to1'][child_haps.index(sample_hap)] += (qc_tbl.at[qchr, 'qARM_START'] - this_aln[2]) / 1000000

                if this_aln[5] >= min_aln_size and separate_single_qry:
                    bases_for_var += write_transmission(out_paf_dir, iden_thresh, alns, qchr)
            else:
                other_bases += this_aln[5]
                samples_bases['Others'][child_haps.index(sample_hap)] += this_aln[5] / 1000000
            contig_sizes.append([qry_tig_len[qchr], sample_hap, '1to1'])
            samples_contig['1to1'][child_haps.index(sample_hap)] += 1
            print(this_paf, file=oneone_writer)
            all_qry_bases += this_aln[5]
            continue

        ## process qry has multiple alignment on ref, which suggests recombination
        is_one2one_contig = True
        valid_paf = []
        for i in range(1, len(sorted_alns)):
            prev_aln, prev_paf = sorted_alns[i - 1]
            curr_aln, curr_paf = sorted_alns[i]
            sorted_out += f'IDEN:{prev_aln[6]},REF:{prev_aln[10]},RSTART:{prev_aln[0]},REND:{prev_aln[1]},QSTART:{prev_aln[2]},QEND:{prev_aln[3]},ORI:{prev_aln[8]};'
            all_qry_bases += prev_aln[5]
            this_contig_aln_bases += prev_aln[5]

            if prev_aln[6] >= iden_thresh and prev_aln[5] >= min_aln_size and prev_aln[4] >= min_aln_size:
                valid_paf.append(sorted_alns[i - 1])

            if i + 1 == len(sorted_alns):
                all_qry_bases += curr_aln[5]
                this_contig_aln_bases += curr_aln[5]
                sorted_out += f'IDEN:{curr_aln[6]},REF:{curr_aln[10]},RSTART:{curr_aln[0]},REND:{curr_aln[1]},QSTART:{curr_aln[2]},QEND:{curr_aln[3]},ORI:{curr_aln[8]};'
                if curr_aln[6] >= iden_thresh and curr_aln[5] >= min_aln_size and prev_aln[4] >= min_aln_size:
                    valid_paf.append(sorted_alns[i])

            ## Check if the same query sequence is aligned to multiple positions
            ro = reciprocal_overlap(prev_aln[2], prev_aln[3], curr_aln[2], curr_aln[3])
            if ro > 0.8:
                is_one2one_contig = False
                if prev_aln[6] >= iden_thresh:
                    onen_bases += prev_aln[5]
                    samples_bases['1toN'][child_haps.index(sample_hap)] += prev_aln[5] / 1000000
                    if i + 1 == len(sorted_alns):
                        other_bases += curr_aln[5]
                        samples_bases['1toN'][child_haps.index(sample_hap)] += curr_aln[5] / 1000000
                else:
                    other_bases += prev_aln[5]
                    samples_bases['Others'][child_haps.index(sample_hap)] += prev_aln[5] / 1000000
                    if i + 1 == len(sorted_alns):
                        other_bases += curr_aln[5]
                        samples_bases['Others'][child_haps.index(sample_hap)] += curr_aln[5] / 1000000

            else:
                if prev_aln[6] >= iden_thresh:
                    oneone_bases += prev_aln[5]
                    samples_bases['1to1'][child_haps.index(sample_hap)] += prev_aln[5] / 1000000
                    if qchr in qc_tbl.index:
                        if prev_aln[3] > qc_tbl.at[qchr, 'qARM_START']:
                            samples_bases_parm['1to1'][child_haps.index(sample_hap)] += (qc_tbl.at[qchr, 'qARM_START'] - prev_aln[2]) / 1000000
                        else:
                            samples_bases_parm['1to1'][child_haps.index(sample_hap)] += prev_aln[5] / 1000000

                    if i + 1 == len(sorted_alns):
                        oneone_bases += curr_aln[5]
                        samples_bases['1to1'][child_haps.index(sample_hap)] += curr_aln[5] / 1000000

                        if qchr in qc_tbl.index:
                            if curr_aln[3] > qc_tbl.at[qchr, 'qARM_START']:
                                samples_bases_parm['1to1'][child_haps.index(sample_hap)] += (qc_tbl.at[qchr, 'qARM_START'] - curr_aln[2]) / 1000000
                            else:
                                samples_bases_parm['1to1'][child_haps.index(sample_hap)] += curr_aln[5] / 1000000


                else:
                    other_bases += prev_aln[5]
                    samples_bases['Others'][child_haps.index(sample_hap)] += prev_aln[5] / 1000000
                    if i + 1 == len(sorted_alns):
                        other_bases += curr_aln[5]
                        samples_bases['Others'][child_haps.index(sample_hap)] += curr_aln[5] / 1000000

        if separate_single_qry:
            bases_for_var += write_transmission(out_paf_dir, iden_thresh, valid_paf, qchr)

        if is_one2one_contig:
            contig_sizes.append([qry_tig_len[qchr], sample_hap, '1to1'])
            samples_contig['1to1'][child_haps.index(sample_hap)] += 1
            for (aln, paf_line) in alns:
                print(paf_line, file=oneone_writer)
        else:
            contig_sizes.append([qry_tig_len[qchr], sample_hap, '1toN'])
            samples_contig['1toN'][child_haps.index(sample_hap)] += 1
            for (aln, paf_line) in alns:
                print(paf_line, file=onen_writer)


        print(f'{qchr}\t{acro_parm_tracker[qchr][0]}\t{len(sorted_alns)}\t{sorted_out[:-1]}', file=aln_writer)

        samples_contig['Total'][child_haps.index(sample_hap)] += 1
        samples_bases['Total'][child_haps.index(sample_hap)] += this_contig_aln_bases/1000000

    print(f'#G3 Contig has 100% identity matches: {perfect_match_num}, bases: {perfect_match_bases/1000000}')

    qry_aligned_bases = round(all_qry_bases/1000000, 2)
    qry_total_oneone_bases = round(oneone_bases / 1000000, 2)
    qry_total_onen_bases = round(onen_bases / 1000000, 2)
    print(f'QRY total aligned bases: {qry_aligned_bases}Mb')
    print(f'QRY total aligned 1to1 bases: {qry_total_oneone_bases}Mb')
    print(f'QRY total aligned 1toN bases: {qry_total_onen_bases}Mb')
    print(f'QRY total aligned Other bases: {round(other_bases/1000000, 2)}Mb')
    print(f'QRY total transmitted bases: {round(bases_for_var / 1000000, 2)}Mb')


    total_qry_contigs = len(qry_alns)
    qry_1to1_contigs = sum(samples_contig['1to1'])
    qry_1toN_contigs = sum(samples_contig['1toN'])

    print(f'Total QRY contigs: {total_qry_contigs}')
    print(f'QRY 1to1 contigs: {qry_1to1_contigs}')
    print(f'QRY 1toN contigs: {qry_1toN_contigs}')
    print(f'QRY Others contigs: {total_qry_contigs - qry_1toN_contigs - qry_1to1_contigs}')

    xticks = np.arange(len(child_haps))
    width = 0.6

    with open(f'{paf_dir}/aligned_bases.json', 'w') as f:
        json.dump(samples_bases, f, indent=4)

    with open(f'{paf_dir}/aligned_bases_parm.json', 'w') as f:
        json.dump(samples_bases_parm, f, indent=4)

    fig, ax = plt.subplots(figsize=(6, 4))
    # other_bases = [i - j - k for i, j, k in zip(samples_bases['Total'], samples_bases['1to1'], samples_bases['1toN'])]
    ax.bar(xticks, samples_bases['Others'], width=width, color='#e08214', label='Others')
    ax.bar(xticks, samples_bases['1to1'], bottom=samples_bases['Others'], width=width, color='#2c78b3', label='1to1')
    ax.bar(xticks, samples_bases['1toN'], bottom=[i + j for i, j in zip(samples_bases['1to1'], samples_bases['Others'])], width=width, color='#32974c', label='1toN')
    # ax.legend(loc='lower left')
    ax.legend('', frameon=False)
    # ax.set_yticks([0, 2000, 4000, 6000])
    # ax.set_yticklabels([0, 2000, 4000, 6000])
    ax.set_ylabel('G3 aligned bases to G2 (Mb)', fontsize=13)
    ax.set_xticks(xticks)
    ax.set_xticklabels(child_haps, rotation=90, fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(f'{paf_dir}/aligned_bases_overview.png', dpi=300)

    # plt.show()


def classify_tranmissions(fai, child_haps, workdir, paf_dir, flag_qc, sd_probe=None, dj_probe=None):

    recombi_tigs = {hap: set() for hap in child_haps}
    total_tigs = {hap: 0 for hap in hap_order}

    for line in open(fai):
        entries = line.strip().split('\t')
        hap = entries[0].split('_')[0] + '_' + entries[0].split('_')[1]
        total_tigs[hap] += 1

    contig_probe = pd.DataFrame()
    contig_dj = pd.DataFrame()
    if sd_probe and dj_probe:
        contig_probe = pd.read_csv(sd_probe, sep='\t', index_col=['Contig'])
        contig_dj = pd.read_csv(dj_probe, sep='\t', index_col=['Contig'])

    contig_qc = {}
    if os.path.exists(f'{flag_qc}/flagger_nucfreq_merged.bed'):
        for line in open(f'{flag_qc}/flagger_nucfreq_merged.bed'):
            entries = line.strip().split('\t')
            if entries[0] not in contig_qc:
                contig_qc[entries[0]] = IntervalTree()
            contig_qc[entries[0]][int(entries[1]): int(entries[2])] = (int(entries[1]), int(entries[2]))

    out_writer = open(f'{workdir}/classified_contig_matches.txt', 'w')
    print('G3-TIG\tG3-FLAGGER\tG3-Probe\tG3-DJ\tG2-TIG\tG2-FLAGGER\tG2-Probe\tG2-DJ\tComb-TAG\tSV-TAG\tSVs', file=out_writer)

    # recomb_view = open(f'{workdir}/recombine_contigs_view.tab', 'w')
    total_matches = 0
    nore, nore_hassv, nore_nosv = 0, 0, 0
    nore_bases, nore_hassv_bases, nore_nosv_bases, reco_bases = 0, 0, 0, 0

    recomb_tracker = {}
    g2_aligned_bases = {}
    for file in os.listdir(paf_dir):
        if 'paf' in file and not file.startswith('.'):
            total_matches += 1
            g3_chr = file.split('.')[0]
            indel_dict = {}
            alned_g2_chrs = set()
            g2_tig_pos, g3_tig_pos = [], []
            g3_aligned_bases = 0
            g2_flag_bases = 0
            g3_flag_bases = 0
            for line in open(f'{paf_dir}/{file}'):
                entries = line.strip().split('\t')
                g3_aligned_bases += (int(entries[3]) - int(entries[2])) / 1000000
                g2_chr = entries[5]

                if g2_chr not in g2_aligned_bases:
                    g2_aligned_bases[g2_chr] = 0
                g2_aligned_bases[g2_chr] += (int(entries[8]) - int(entries[7])) / 1000000

                g2_flag_ovlps = contig_qc[entries[5]].overlap(int(entries[7]), int(entries[8])) if entries[5] in contig_qc else []
                if g2_flag_ovlps:
                    for ovlp in g2_flag_ovlps:
                        g2_flag_bases += ovlp.end - ovlp.begin

                g3_flag_ovlps = contig_qc[entries[0]].overlap(int(entries[2]), int(entries[3])) if entries[0] in contig_qc else []
                if g3_flag_ovlps:
                    for ovlp in g3_flag_ovlps:
                        g3_flag_bases += ovlp.end - ovlp.begin

                g2_tig_pos.append(f'{entries[5]}:{entries[7]}-{entries[8]}')
                g3_tig_pos.append(f'{entries[0]}:{entries[2]}-{entries[3]}')

                if g2_chr not in indel_dict:
                    indel_dict[g2_chr] = []

                alned_g2_chrs.add(g2_chr)
                cg = [i.split(":")[-1] for i in entries[12:] if i[:2] == 'cg'][0]

                cg_tuple = cigar_tuple(cg)
                indels = insdel_from_cigar(cg_tuple)
                for (ref_pos, read_pos, size, svtype) in indels:
                    if size < 50:
                        continue
                    var_id = f'{ref_pos + int(entries[7])}-{size}-{svtype}'
                    indel_dict[g2_chr].append(var_id)
            var_list = []
            tag_1 = 'Non-Recomb' if len(alned_g2_chrs) == 1 else 'Recomb'
            tag_2 = 'No-SVs'
            out_g2_chr = []
            out_g2_chr_probe = []
            out_g2_chr_dj = []


            for g2_chr, vars in indel_dict.items():
                this_vars = f'{g2_chr}:'
                out_g2_chr.append(g2_chr)
                # out_g2_chr_flag.append(str(contig_qc[g2_chr]) if len(contig_qc)!=0 else '.')

                if g2_chr in contig_probe.index:
                    hits = contig_probe.at[g2_chr, 'Probe']
                    if type(hits) == str:
                        out_g2_chr_probe.append(hits)
                    else:
                        out_g2_chr_probe.append('-'.join(hits.tolist()))
                else:
                    out_g2_chr_probe.append('.')
                if g2_chr in contig_dj.index:
                    hits = contig_dj.at[g2_chr, 'Probe']
                    if type(hits) == str:
                        out_g2_chr_dj.append(hits)
                    else:
                        out_g2_chr_dj.append('-'.join(hits.tolist()))
                else:
                    out_g2_chr_dj.append('.')

                if len(vars) > 0:
                    tag_2 = 'Has-SVs'
                    for var in vars:
                        this_vars += f'{var},'
                    var_list.append(this_vars[:-1])
            if tag_1 == 'Non-Recomb':
                nore += 1
                nore_bases += g3_aligned_bases
                if tag_2 == 'Has-SVs':
                    nore_hassv += 1
                    nore_hassv_bases += g3_aligned_bases
                else:
                    nore_nosv += 1
                    nore_nosv_bases += g3_aligned_bases
            else:
                ## debug usage to check recombine
                # png_name = f'{g3_chr}.iden99.aln1Mb.png'
                # print(f'{file}\t{png_name}', file=recomb_view)

                reco_bases += g3_aligned_bases
                g3_hap = g3_chr.split('_')[0] + '_' + g3_chr.split('_')[1]
                recombi_tigs[g3_hap].add(g3_chr)
                froze_g2_chr = frozenset(out_g2_chr)
                if froze_g2_chr not in recomb_tracker:
                    recomb_tracker[froze_g2_chr] = []
                recomb_tracker[froze_g2_chr].append([g3_chr, g3_aligned_bases])

            ## Add percentage of bad bases for each contig
            # bad_pcrt.append([contig_qc[g3_chr] * 100, g3_chr.split('-')[0], 'G3'])
            # bad_pcrt.append([contig_qc[out_g2_chr] * 100, out_g2_chr.split('-')[0], 'G2'])

            g2_tigs = ';'.join(out_g2_chr)
            # g2_flag = ';'.join(out_g2_chr_flag)
            g2_probe = ';'.join(out_g2_chr_probe)
            g2_dj = ';'.join(out_g2_chr_dj)

            g3_probe, g3_dj = '.', '.'
            if g3_chr in contig_probe.index:
                hits = contig_probe.at[g3_chr, 'Probe']
                if type(hits) == str:
                    g3_probe = hits
                else:
                    g3_probe = '-'.join(hits.tolist())

            if g3_chr in contig_dj.index:
                hits = contig_dj.at[g3_chr, 'Probe']
                if type(hits) == str:
                    g3_dj = hits
                else:
                    g3_dj = '-'.join(hits.tolist())

            # g3_contig_qc = contig_qc[g3_chr] if len(contig_qc) > 0 else '.'
            var_out = '.' if len(var_list) == 0 else ','.join(var_list)
            print(f'{g3_chr}\t{g3_flag_bases}\t{g3_probe}\t{g3_dj}\t{g2_tigs}\t{g2_flag_bases}\t{g2_probe}\t{g2_dj}\t{tag_1}\t{tag_2}\t{var_out}', file=out_writer)

    track_out = open(f'{workdir}/recomb_tracker.tsv', 'w')
    g2_tig_counts = {}
    print('G2\tG3-Num\tG3-Bases\tG3', file=track_out)
    for g2_chr, g3_list in recomb_tracker.items():
        if len(list(g2_chr)) not in g2_tig_counts:
            g2_tig_counts[len(list(g2_chr))] = 0
        g2_tig_counts[len(list(g2_chr))] += 1

        g3_tigs = [ele[0] for ele in g3_list]

        g3_bases = sum([ele[1] for ele in g3_list])
        g2_str = ','.join(list(g2_chr))
        g3_str = ','.join(g3_tigs)
        print(f'{g2_str}\t{len(g3_tigs)}\t{g3_bases}\t{g3_str}', file=track_out)

    print("G3 bases: ",reco_bases, nore_bases)
    print("G2 bases: ", sum(g2_aligned_bases.values()))

def write_transmission(out_paf_dir, iden_thresh, valid_paf_lines, qchr):
    transmitted_bases = 0
    if len(valid_paf_lines) > 0:
        qry_paf = open(f'{out_paf_dir}/qry_contig_pafs/{qchr}.iden{iden_thresh}.aln1Mb.paf', 'w')
        qry_tsv = open(f'{out_paf_dir}/qry_contig_pafs/{qchr}.iden{iden_thresh}.aln1Mb.tsv', 'w')
        print('qry\tqry_start\tqry_end\tref\tref_start\tref_end', file=qry_tsv)
        for (aln, paf_line) in valid_paf_lines:
            print(paf_line, file=qry_paf)
            entries = paf_line.strip().split('\t')
            transmitted_bases += int(entries[3]) - int(entries[2])
            print(f'{entries[0]}\t{entries[2]}\t{entries[3]}\t{entries[5]}\t{entries[7]}\t{entries[8]}', file=qry_tsv)

        qry_paf.close()
        qry_tsv.close()

    return transmitted_bases


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_autopct

def reciprocal_overlap(begin_a, end_a, begin_b, end_b):

    overlap = min(end_a, end_b) - max(begin_a, begin_b)

    return round(min([overlap / (end_a - begin_a), overlap / (end_b - begin_b)]), 3)
