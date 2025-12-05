#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/4/25
'''
import os

from ParsePaf import classify_aligns, classify_tranmissions


def main():
    g3_child = ['NA12879_1', 'NA12879_2', 'NA12881_1', 'NA12881_2',
                         'NA12882_1', 'NA12882_2', 'NA12883_1', 'NA12883_2',
                         'NA12884_1', 'NA12884_2', 'NA12885_1', 'NA12885_2',
                         'NA12886_1', 'NA12886_2', 'NA12887_1', 'NA12887_2']

    g4_fam1_child = ['K200081_1', 'K200081_2', 'K200082_1', 'K200082_2', 'K200084_1', 'K200084_2', 'K200085_1',
                        'K200085_2', 'K200086_1', 'K200086_2', 'K200087_1', 'K200087_2']

    wfmash_dir = './'
    contig_qc_dir = './annotation/'

    trans_dir = f'{wfmash_dir}/classify_alns'
    if not os.path.exists(trans_dir):
        os.mkdir(trans_dir)
        os.mkdir(f'{trans_dir}/qry_contig_pafs')


    ########################################################################
    #### Step1: Classify offspring to parent alignment to 1to1 and 1toN ####
    ########################################################################
    classify_aligns(f'{wfmash_dir}/sequence.size', g4_fam1_child, wfmash_dir, f'{wfmash_dir}/g2g3.paf',
                       f'{wfmash_dir}/g2g3.qry_aligns.txt', trans_dir,
                         99, 1000000, separate_single_qry=True)

    ###############################################################################
    #### Step2: Find transmission and recombination from classified alignments ####
    ###############################################################################

    classify_tranmissions(f'{wfmash_dir}/sequence.size', g4_fam1_child, f'{wfmash_dir}/classify_alns',
        f'{trans_dir}/qry_contig_pafs', contig_qc_dir)

if __name__ == '__main__':
    main()
