#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/6/25
'''
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from scipy import stats
from statsmodels.stats.multitest import multipletests

all_children = ['NA12879', 'NA12881', 'NA12882', 'NA12883', 'NA12884','NA12885', 'NA12886', 'NA12887', 'K200081', 'K200082', 'K200084', 'K200085', 'K200086', 'K200087']
SAT_COLORS = {'bSat': '#FFD63A', 'HSat1A': '#00DE60', 'HSat1B': '#1B998B','HSat3': '#fdb462', 'ACRO': '#3E7C17', 'aSat': '#990000', 'PHR_LOC': '#FF6600',
              'Unlabelled': '#bdbdbd', 'Others': '#BED7DC', 'DJ': '#990099', 'DJflank': '#990099', 'DJarm': '#990099', 'PJ': '#FF6666',
              'SegDup': '#386cb0', 'rDNA': '#9966FF', 'SST1': '#FF2DD1', 'parm-sd': '#386cb0', 'qarm-sd': '#386cb0', 'SEQ': '#f0f0f0',
              'SD': '#386cb0', 'aSat_HOR': '#9970ab', 'CER': '#00CCCC', 'GAP': '#252525'}

def figure6a():

    rate_list = []
    chrom_rate_list = []

    dnm_tbl = pd.read_csv('./final_denovo_SNPs.tsv', sep='\t', header=0)
    dnm_tbl['mutation_count'] = [1 for _ in range(len(dnm_tbl))]

    dnm_tbl['hap'] = dnm_tbl.apply(lambda row: row['child'] + '_Pat' if row['parent_type'] == 'father' else row['child'] + '_Mat', axis=1)
    dnm_sum_byannot_sample = dnm_tbl.groupby(['repeats', 'hap'])['mutation_count'].sum().reset_index()
    dnm_sum_byannot_sample.loc[len(dnm_sum_byannot_sample)] = {'repeats': 'SEQ', 'hap': 'NA12879_Mat', 'mutation_count': 0}
    dnm_sum_byannot_sample.loc[len(dnm_sum_byannot_sample)] = {'repeats': 'SEQ', 'hap': 'NA12881_Pat','mutation_count': 0}

    dnm_sum_by_sample = dnm_tbl.groupby(['hap'])['mutation_count'].sum().reset_index()


    ## Get callable bases
    subset_dict = {'wfmash_g2g3': 'G2-G3', 'wfmash_g3_NA12879': 'G4_fam1'}
    for subset, fam in subset_dict.items():
        df_rate = pd.read_csv(f'{subset}_final_dnms_rate.tsv', sep='\t', header=[0])
        df_rate['bases_mb'] = df_rate.apply(lambda row: row['bases'] / 1000000, axis=1)
        df_rate['fam'] = [fam for _ in range(len(df_rate))]

        rate_list.append(df_rate)

        df_chrom = pd.read_csv(f'{subset}_final_dnms_rate_bychrom.tsv', sep='\t', header=[0])
        df_chrom['fam'] = [fam for _ in range(len(df_chrom))]
        chrom_rate_list.append(df_chrom)

    all_chrom_rate = pd.concat(chrom_rate_list)
    all_rate = pd.concat(rate_list)
    all_rate['hap_index'] = all_rate.apply(lambda row: '{0}_{1}'.format(row['g3_sample'], row['hap']), axis=1)
    all_rate.set_index('hap_index', inplace=True)
    dnm_sum_by_sample['callable'] = dnm_sum_by_sample.apply(lambda row: all_rate.at[row['hap'], 'bases'] if row['hap'] in all_rate.index else 0, axis=1)
    dnm_sum_by_sample['manual_rate'] = dnm_sum_by_sample.apply(lambda row: row['mutation_count'] / row['callable'], axis=1)

    dnm_sum_by_sample.loc[len(dnm_sum_by_sample)] = {'callable': 0, 'hap': 'NA12879_Mat', 'mutation_count': 0, 'manual_rate': 0}
    dnm_sum_by_sample.loc[len(dnm_sum_by_sample)] = {'callable': 0, 'hap': 'NA12881_Pat','mutation_count': 0, 'manual_rate': 0}

    mat_dnm = dnm_sum_by_sample.loc[dnm_sum_by_sample['hap'].str.contains('Mat')]
    pat_dnm = dnm_sum_by_sample.loc[dnm_sum_by_sample['hap'].str.contains('Pat')]

    # dnm_sum_by_sample.to_csv()

    mat_avg_rate = mat_dnm['manual_rate'].mean()
    pat_avg_rate = pat_dnm['manual_rate'].mean()
    print('All mutation rate: ', all_rate['manual_rate'].mean())
    print('Maternal mutation rate: ', mat_dnm['manual_rate'].mean())
    print('Paternal mutation rate: ', pat_dnm['manual_rate'].mean())

    # chrom_rate_mean = all_chrom_rate.groupby(['chrom', 'hap'])['manual_rate'].mean().reset_index()

    annot_type = dnm_sum_byannot_sample['repeats'].unique().tolist()

    parent_order = ['Pat', 'Mat']
    sorted_hap_order = sorted(dnm_sum_byannot_sample['hap'].unique().tolist(), key=lambda x: (parent_order.index(x.split('_')[1]), all_children.index(x.split('_')[0])))

    sorted_mut_rate = []
    sorted_pat_rate = []
    sorted_mat_rate = []
    all_rate.reset_index(drop=True, inplace=True)
    for i, hap in enumerate(sorted_hap_order):
        qry_sample, tmp_hap = hap.split('_')
        mut_rate = dnm_sum_by_sample.loc[dnm_sum_by_sample['hap']==hap]['manual_rate'].values[0]
        if tmp_hap == 'Pat':
            sorted_pat_rate.append(mut_rate)
        elif tmp_hap == 'Mat':
            sorted_mat_rate.append(mut_rate)
        sorted_mut_rate.append(mut_rate)

    legends = []

    width=0.6

    pat_dnm = dnm_sum_byannot_sample.loc[dnm_sum_byannot_sample['hap'].str.contains('Pat')]
    mat_dnm = dnm_sum_byannot_sample.loc[dnm_sum_byannot_sample['hap'].str.contains('Mat')]
    fig, axes = plt.subplots(4,1, sharex=True, figsize=(6, 8))

    xticks = np.arange(len(all_children) - 1)
    for j, dnm_counts in enumerate([pat_dnm, mat_dnm]):
        bottom_values = []
        sample_order = sorted(dnm_counts['hap'].unique().tolist(), key=lambda x: all_children.index(x.split('_')[0]))
        for i, sat in enumerate(annot_type):
            this_sat_num = dnm_counts.loc[dnm_counts['repeats'] == sat]
            print(f'{sat}, #SNV:', this_sat_num['mutation_count'].sum())
            curr_vals = [0 for _ in range(len(xticks))]

            for idx, row in this_sat_num.iterrows():
                hap_index = sample_order.index(row['hap'])
                curr_vals[hap_index] += row['mutation_count']

            legends.append(Patch(label=sat, color=SAT_COLORS[sat]))

            if i == 0:
                axes[j].bar(xticks, curr_vals, color=SAT_COLORS[sat], width=width)
            else:
                current_bottoms=[sum(elements) for elements in zip(*bottom_values[0:i])]
                axes[j].bar(xticks, curr_vals, bottom=current_bottoms, width=width, color=SAT_COLORS[sat])

            bottom_values.append(curr_vals)
        axes[j].set_ylabel('DNM counts', fontsize=14)
        axes[j].set_yticks([0, 3, 6, 9], labels=[0, 3, 6, 9], fontsize=13)

    axes[2].plot(xticks, sorted_pat_rate, lw=2, marker='o', color='#3288bd')
    axes[2].plot(xticks, [pat_avg_rate for _ in range(len(xticks))], lw=1.5, ls='--',color='#3288bd')
    axes[3].plot(xticks, sorted_mat_rate, lw=2, marker='o', color='#4d9221')
    axes[3].plot(xticks, [mat_avg_rate for _ in range(len(xticks))], lw=1.5, ls='--', color='#4d9221')
    axes[3].set_xticks(xticks, labels=['NA12879', 'NA12881', 'NA12882', 'NA12883', 'NA12884','NA12885', 'NA12886', 'NA12887', '200081', '200082', '200084', '200085', '200087'], rotation=90, fontsize=14)

    for k in [2, 3]:
        axes[k].set_ylabel('DNM rate (10^-7)', fontsize=14)
        axes[k].set_yticks([0.5e-7, 1.5e-7, 2.5e-7], [0.5, 1.5, 2.5] , fontsize=13)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    plt.show()


def update_chrom_mut_rate(row, df):
    return row['mutation_count']/df.loc[row['sample'], row['region']]['callable_space']
def update_cen_mut_rate(row, df):
    if row['region'] == 'total_cen':
        return row['mutation_count']/df.at[row['sample'], 'total']
    elif row['region'] == 'HOR':
        return row['mutation_count'] / df.at[row['sample'], 'HOR']
    elif row['region'] == 'flanking':
        return row['mutation_count'] / df.at[row['sample'], 'flanking']

def figure6c():
    rate_to_plot = []
    ## autosome variant
    auto_space = pd.read_csv('./denovo/mutation_rates/callable_regions.tsv', sep='\t', index_col=['sample', 'region'])
    auto_counts = pd.read_csv('./denovo/mutation_rates/dnm_counts/final_dnm_pzm_counts.tsv', sep='\t', header=0)
    auto_counts['fam'] = ['all' for _ in range(len(auto_counts))]
    auto_counts['rate'] = auto_counts.apply(lambda row: update_chrom_mut_rate(row, auto_space), axis=1)
    rate_to_plot.append(auto_counts.loc[(auto_counts['region']=='autosomes')|(auto_counts['region']=='all_segdup')][['region', 'rate', 'fam']])

    ## chrY mutation
    chrY_rate = pd.read_csv('./denovo/mutation_rates/chrY.tsv', sep='\t', header=0)
    chrY_rate.rename(columns={'mutation_rate': 'rate'}, inplace=True)
    chrY_rate['fam'] = ['all' for _ in range(len(chrY_rate))]
    rate_to_plot.append(chrY_rate[['region', 'rate', 'fam']])

    ## cen mutation
    cen_counts = pd.read_csv('./denovo/mutation_rates/centromeres/dnm_counts_by_sample.tsv', sep='\t', header=0)
    cen_space = pd.read_csv('./denovo/mutation_rates/centromeres/callable_regions.tsv', sep='\t', index_col=['sample'])

    cen_counts['rate'] = cen_counts.apply(lambda row: update_cen_mut_rate(row, cen_space), axis=1)
    cen_counts['fam'] = ['all' for _ in range(len(cen_counts))]
    rate_to_plot.append(cen_counts.loc[cen_counts['region']=='total_cen'][['region', 'rate', 'fam']])

    ## acro final denovo SNPs
    dnm_tbl = pd.read_csv(f'./final_denovo_SNPs.tsv', sep='\t', header=0)
    dnm_tbl['mutation_count'] = [1 for _ in range(len(dnm_tbl))]
    dnm_sum_by_sample = dnm_tbl.groupby(['child'])['mutation_count'].sum().reset_index()

    ## acro mutation rate
    acro_rates = []
    subset_dict = {'wfmash_g2g3': 'G2-G3', 'wfmash_g3_NA12879': 'G4_fam1'}
    for subset, fam in subset_dict.items():
        acro_rate_all = pd.read_csv(f'{subset}_final_dnms_rate.tsv', sep='\t', header=0)
        for sample, group in acro_rate_all.groupby('g3_sample'):
            # mut_num = group['manual_count'].sum()
            mut_num = dnm_sum_by_sample.loc[dnm_sum_by_sample['child']==sample]['mutation_count'].values[0]
            mut_space = group['bases'].sum()
            rate = 0 if mut_num == 0 else mut_num/mut_space
            acro_rates.append(['all_acro', rate, fam])
    df_acro_rates = pd.DataFrame(acro_rates, columns=['region', 'rate', 'fam'])

    rate_to_plot.append(df_acro_rates)

    all_rates = pd.concat(rate_to_plot)
    all_rates.reset_index(drop=True, inplace=True)

    # fig, axes = plt.subplots(1,2, sharey=True, figsize=(9, 4))
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.stripplot(data=all_rates, x='region', y='rate', hue='fam', order=['autosomes', 'all_segdup', 'total_cen', 'chrY', 'all_acro'], size=7, ax=ax)
    sns.boxplot(x="region", y="rate", data=all_rates, color='#ffffff', order=['autosomes', 'all_segdup', 'total_cen', 'chrY', 'all_acro'], ax=ax)
    # sns.stripplot(data=all_rates, x='region', y='rate', hue='fam', order=['chr13', 'chr14', 'chr15', 'chr21', 'chr22'], size=7, ax=axes[1])
    # sns.boxplot(x="region", y="rate", data=all_rates, color='#ffffff',order=['chr13', 'chr14', 'chr15', 'chr21', 'chr22'], ax=axes[1])
    ax.set_ylabel('Mutation rate (10^-7)', fontsize=15)
    ax.set_yticks([0, 1e-7, 2e-7, 3e-7, 4e-7], labels=[0, 1, 2, 3, 4], fontsize=13)
    ax.set_xticks(np.arange(5), labels=['Auto', 'SegDup', 'CEN', 'Yq12', 'Acro'], fontsize=14)

    ax.set_xlabel('')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    plt.show()


def figure6d():
    ## acro final denovo SNPs
    dnm_tbl = pd.read_csv('./final_denovo_SNPs.tsv', sep='\t', header=0)
    dnm_tbl.set_index('id', inplace=True)
    ## original platinum call spectrum
    autosome_calls = pd.read_csv('./denovo/final_calls/final_calls_without_false_pzms/final_autosomal_snvs.tsv', sep='\t', header=0)
    autosome_calls = autosome_calls.loc[autosome_calls['variant_type']!='pzm']
    autosome_calls['region'] = ['autosome' for _ in range(len(autosome_calls))]

    ## acrocentric call spectrum
    acro_call_list = []
    acro_spec_tbl = pd.read_csv('./spectrum_var.tab', sep='\t', header=0)
    for idx, row in acro_spec_tbl.iterrows():
        sample, parent_type, region = row['sample'], row['parent_type'], row['region']
        this_spec = pd.read_csv(f'./spectrum/annotated_snvs/{sample}_from_{parent_type}.{region}.snvs.tsv', sep='\t', header=0)
        this_spec = this_spec.loc[this_spec['id'].isin(dnm_tbl.index)]
        if len(this_spec) > 0:
            this_spec['region'] = ['acrocentric' for _ in range(len(this_spec))]
            acro_call_list.append(this_spec)
    acro_calls = pd.concat(acro_call_list)

    all_mut = pd.concat([autosome_calls[['region', 'triplet', 'CpG', 'substitution_type']], acro_calls[['region', 'triplet', 'CpG', 'substitution_type']]])
    all_mut['five_prime'] = all_mut.apply(lambda row: row['triplet'][0], axis=1)
    all_mut['three_prime'] = all_mut.apply(lambda row: row['triplet'][2], axis=1)
    all_mut['mutation'] = all_mut.apply(lambda row: '{0}>{1}'.format(row['triplet'][1], row['triplet'][4]), axis=1)

    mutation_spectrum = all_mut.groupby(['region', 'mutation', 'three_prime', 'substitution_type']).count().reset_index()[['region', 'mutation', 'three_prime', 'substitution_type','triplet']]
    mutation_spectrum.rename(columns={'triplet': 'mutation_count'}, inplace=True)
    mutation_spectrum_region_sum = mutation_spectrum.groupby(['region'])['mutation_count'].sum().reset_index()
    mutation_spectrum_region_sum.set_index('region', inplace=True)

    mutation_spectrum['fraction'] = mutation_spectrum.apply(lambda row: row['mutation_count']/mutation_spectrum_region_sum.at[row['region'], 'mutation_count'], axis=1)
    dinu_spec = []
    for (region, mutation, three_prime), group in mutation_spectrum.groupby(['region', 'mutation', 'three_prime']):
        mut_frac = group['fraction'].sum()
        mut_count = group['mutation_count'].sum()
        dinu_spec.append([region, mutation, three_prime, mut_frac, mut_count])

    df_dinu_spec = pd.DataFrame(dinu_spec, columns=['region', 'mutation', 'three_prime', 'fraction', 'count'])
    mutation_types = df_dinu_spec['mutation'].unique().tolist()

    fig, axes = plt.subplots(1, len(mutation_types), sharey=True, sharex=True, figsize=(11, 4))

    for i, mut_type in enumerate(mutation_types):
        ax = axes[i]
        sub_spec = df_dinu_spec.loc[df_dinu_spec['mutation'] == mut_type]
        sns.barplot(data=sub_spec, x='three_prime', y='fraction', hue='region', hue_order=['autosome', 'acrocentric'], order=['A', 'C', 'G', 'T'], ax=ax)
        mut_a, mut_b = mut_type.split('>')

        ax.set_title(f'{mut_a}pN>{mut_b}pN', fontsize=13)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_xticks([0, 1, 2, 3], labels=['A', 'C', 'G', 'T'], fontsize=13)
        ax.set_xlabel('')
        ax.legend('', frameon=False)
    axes[0].set_ylabel('Fraction of DNM calls', fontsize=15)
    fig.tight_layout()
    plt.show()

    ## test the significance of mutations
    total_mutations = df_dinu_spec.groupby(['region'])['count'].sum().reset_index()
    total_mutations.rename(columns={'count': 'total'}, inplace=True)

    mut_base_tbl = df_dinu_spec.groupby(['mutation', 'region'])['count'].sum().reset_index()
    mut_types = mut_base_tbl['mutation'].unique().tolist()
    pvals = []

    for mut in mut_types:
        this_mut = mut_base_tbl.loc[mut_base_tbl['mutation'] == mut]
        test_tbl = pd.merge(this_mut, total_mutations, on='region', how='left')
        test_tbl['total'] = test_tbl.apply(lambda row: row['total'] - row['count'], axis=1)
        stats, p = fisher_exact(test_tbl[['count', 'total']].values.tolist(), alternative='two-sided')

        pvals.append([mut, p])
        if p < 0.05:
            # print('-----{0} {1}-----'.format(mut.split('_')[0], mut.split('_')[1]))
            print(f'{mut}\t{p}')

    sorted_pvals = sorted(pvals, key=lambda x:x[1])
    df_pvals = pd.DataFrame(sorted_pvals, columns=['mut', 'pval'])
    df_pvals.to_csv('spectrum/fisher_mut_test.tsv', sep='\t', header=True, index=False)
    qvals = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
    print(qvals)

    ## test the significance of dinucleotide mutations
    mutation_tbl = df_dinu_spec.copy()
    mutation_tbl['mutation'] = df_dinu_spec.apply(lambda row: '{0}_{1}'.format(row['mutation'], row['three_prime']),                                                  axis=1)
    mut_types = mutation_tbl['mutation'].unique().tolist()
    pvals = []
    total_mutations = mutation_tbl.groupby(['region'])['count'].sum().reset_index()
    total_mutations.rename(columns={'count': 'total'}, inplace=True)
    for mut in mut_types:
        this_mut = mutation_tbl.loc[mutation_tbl['mutation'] == mut]
        test_tbl = pd.merge(this_mut, total_mutations, on='region', how='left')
        test_tbl['total'] = test_tbl.apply(lambda row: row['total'] - row['count'], axis=1)

        stats, p = fisher_exact(test_tbl[['count', 'total']].values.tolist(), alternative='two-sided')

        pvals.append(p)
        if p < 0.05:
            print('-----{0} {1}-----'.format(mut.split('_')[0], mut.split('_')[1]))
            print(f'\t{p}')

    qvals = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
    print(qvals)



def main():

    figure6a()
    figure6c()
    figure6d()

if __name__ == '__main__':
    main()
