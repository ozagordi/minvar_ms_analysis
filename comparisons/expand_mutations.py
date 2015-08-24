#!/usr/bin/env python3.4
'''Adds the wild type frequency to the mutations detected with MinVar or
VirVarSeq.'''
import sys
import pandas as pd
from Bio import SeqIO

muts = pd.read_csv(sys.argv[1])

# if VirVarSeq, change names
if 'POSITION' in muts.columns:
    muts.drop(['REF_AA','POSITION'], 1, inplace=True)
    muts.rename(columns={'AA': 'mut', 'FREQ': 'freq'}, inplace=True)
    muts.replace(to_replace='.', value='-', inplace=True)
# if true mixes
if 'wt' in muts.columns:
    muts.drop('wt', 1, inplace=True)
muts = muts[muts.gene.isin(["RT", "protease"])]
protease = list(SeqIO.parse('protease.faa', 'fasta'))[0].seq
RT = list(SeqIO.parse('RT.faa', 'fasta'))[0].seq

assert str(protease[53]) == 'I', protease[53]
assert str(RT[58]) == 'P'

f = muts.groupby(['gene', 'pos'])
for gp, group in f:
    gene, pos = gp
    mf = group.sum().freq
    if abs(1. - mf) < 1.E-3 or gene in ['GagPolTF', 'integrase']:
        continue
    if gene == 'protease':
        wt = protease[pos - 1]
    elif gene == 'RT':
        wt = RT[pos - 1]
    fwt = 1.0 - mf
    dict_here = {'gene': gene, 'pos': pos, 'mut': str(wt), 'freq': fwt}
    muts = muts.append(dict_here, ignore_index=True)
muts = muts.sort(columns=['gene', 'pos'])
muts.to_csv(sys.stdout, sep=',', float_format='%6.4f', index=False)
