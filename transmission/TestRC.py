#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/6/25
'''

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

CHROMS = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']

def test_rc_bp(test_range):
    print(test_range)


    chrom_size = pd.read_csv(f'./rc_test/qarm_vs_qarm/chm13_qarm.bed',sep='\t', names=['chrom', 'start', 'end'])
    chrom_size.set_index('chrom', inplace=True)
    ## This test randomly samples 1000 5Mb regions on the rest of the autosome regions to build the distribution
    ## excluded one chr21 parm ectopic recombination event
    acro_counts = {'chr13': 4, 'chr15': 3, 'chr14': 6, 'chr21': 3, 'chr22': 2}
    perm_cnt = []
    perm_cnt_chrom = []
    for i in range(1, 5001):
        shl_idx = pd.read_csv(f'./rc_test/{test_range}/tmp/shuffle_{i}_bp_cnt.txt', sep='\t',
                              names=['chrom', 'start', 'end', 'count'])
        perm_cnt.append([i, shl_idx['count'].sum()])
        perm_cnt_chrom.append(shl_idx)
    df_perm_cnt = pd.DataFrame(perm_cnt, columns=['perm_idx', 'count'])

    tmp_perm_chrom = pd.concat(perm_cnt_chrom)
    legend = [Line2D([0], [0], color='r', ls='--', label='Observation')]

    fig, ax = plt.subplots(figsize=(7, 4))
    # sns.kdeplot(data=df_perm_cnt, x='count', ax=ax)
    sns.histplot(data=df_perm_cnt, x='count', binwidth=1, kde=True, ax=ax)
    ax.set_xlabel('#Recombination events', fontsize=14)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylabel('#Randomly shuffled regions', fontsize=14)
    if test_range == 'qarm_vs_all':
        ax.axvline(x=18, lw=1.5, ls='--', color='red')
        ax.set_title('Euchromatin null distribution', fontsize=14)
    elif test_range == 'parm_vs_parm':
        ax.axvline(x=0, lw=1.5, ls='--', color='red')
        ax.set_title('parm null distribution', fontsize=14)
    fig.tight_layout()
    ax.legend(handles=legend)

    data = df_perm_cnt['count'].to_numpy()
    if test_range == 'qarm_vs_all':

        empirical_p_value = np.mean(data >= 18)
        print(empirical_p_value, df_perm_cnt['count'].mean(), df_perm_cnt['count'].median())

        for chrom in CHROMS:
            chrom_cnt = tmp_perm_chrom.loc[tmp_perm_chrom['chrom']==chrom]['count'].to_numpy()
            chrom_empirical_p = np.mean(chrom_cnt >= acro_counts[chrom])
            print(f'{chrom}: {chrom_empirical_p}', acro_counts[chrom], tmp_perm_chrom.loc[tmp_perm_chrom['chrom']==chrom]['count'].mean())

    elif test_range == 'parm_vs_parm':
        empirical_p_value = np.mean(data ==0)
        print(empirical_p_value, df_perm_cnt['count'].mean())

    plt.show()

def main():

    test_rc_bp('qarm_vs_all')
    test_rc_bp('parm_vs_parm')


if __name__ == '__main__':
    main()
