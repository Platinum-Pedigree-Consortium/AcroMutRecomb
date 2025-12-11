#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/6/25
'''
import pandas as pd
from pyfaidx import Fasta
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import random

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def normalize_pyrimidine(ref, alt, tri):
    ref, alt, tri = ref.upper(), alt.upper(), tri.upper()
    if ref in ("A","G"):
        return revcomp(ref), revcomp(alt), revcomp(tri)
    return ref, alt, tri

def compute_mut_contexts(fasta, muts, k):
    half = k // 2
    nctxs = []
    nsubs = []
    for i, row in muts.iterrows():
        chrom = row['chrom']
        pos = int(row['pos'])
        try:
            seq = fasta[chrom][pos - half - 1:pos + half].seq.upper()
        except Exception:
            nctxs.append(None)
            nsubs.append(None)
            continue
        if len(seq) != k or 'N' in seq:
            nctxs.append(None)
            nsubs.append(None)
            continue
        nref, nalt, nctx = normalize_pyrimidine(row['ref'], row['alt'], seq)
        nctxs.append(nctx)
        nsubs.append(f"{nref}>{nalt}")
    muts['nctx'] = nctxs
    muts['nsub'] = nsubs
    muts = muts.dropna(subset=['nctx', 'nsub']).copy()
    return muts


def load_mutations(fn):
    df = pd.read_csv(fn, sep='\t', header=0, usecols=[0,1,2,3], names=['chrom','pos','ref','alt'], dtype={'chrom':str})
    return df

def main():
    workdir = './dnm_all'
    ## mutation to test

    sub = 'A>C'
    muts = pd.read_csv(f'{workdir}/final_denovo_SNPs.bed', sep='\t', header=0, names=['chrom','pos','ref','alt'])
    muts = compute_mut_contexts(Fasta(f'{workdir}/final_denovo_REF.fa'), muts, 3)
    sub_ref, sub_alt = sub.split(">")
    if sub_ref in ("A", "G"):
        # map to pyrimidine ref/alt by reverse complement of single base (no context needed)
        sub_ref, sub_alt = revcomp(sub_ref), revcomp(sub_alt)

    requested_sub = f"{sub_ref}>{sub_alt}"
    muts_sub = muts[muts['nsub'] == requested_sub].copy()
    N = len(muts_sub)

    obs_counter = Counter(muts_sub['nctx'].tolist())

    ## Create background context
    bck_muts = pd.read_csv(f'{workdir}/michelle_dnms/fully_annotated_autosomal_snvs.bed', sep='\t', usecols=[0,1,3,4,5], names=['chrom','pos','ref','alt', 'spe'])
    bg_counts = Counter()
    for idx, row in bck_muts.iterrows():
        tri = row['spe'].split('-')[0]
        center = tri[1]
        if center in ('A', 'G'):
            tri = revcomp(tri)
        bg_counts[tri] += 1

    B = sum(bg_counts.values())
    contexts_to_test = sorted([k for k in bg_counts.keys()])
    for ctx in contexts_to_test:
        n_i = obs_counter.get(ctx, 0)
        b_i = bg_counts.get(ctx, 0)
        lam = N * (b_i / B) if B > 0 else 0.0
        # one-sided test (enrichment)
        p = poisson.sf(n_i - 1, lam) if lam > 0 else (1.0 if n_i == 0 else 0.0)
        fold = (n_i / lam) if lam > 0 else (np.nan if n_i == 0 else np.inf)
        if p < 0.05:
            print(ctx, p, fold)

if __name__ == '__main__':
    main()
