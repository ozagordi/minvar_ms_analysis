#!/usr/bin/env python3.4
import sys
import csv
import pandas as pd
from Bio import SeqIO

# Sanger results
sanger_file = 'Minority24.tsv'
with open(sanger_file) as csvfile:
    pat_list = csv.reader(csvfile, delimiter='\t')
    for p in pat_list:
        patient, sm_prot, sm_RT, subtype = p
        if patient == sys.argv[1]:
            sanger_muts_protease = {}
            sanger_muts_RT = {}
            for m in sm_prot.split(','):
                try:
                    pos, mut = int(m[:-1]), m[-1]
                except ValueError:
                    continue
                try:
                    sanger_muts_protease[pos].append(mut)
                except KeyError:
                    sanger_muts_protease[pos] = [mut]
            for m in sm_RT.split(','):
                pos, mut = int(m[:-1]), m[-1]
                try:
                    sanger_muts_RT[pos].append(mut)
                except KeyError:
                    sanger_muts_RT[pos] = [mut]
            break

sanger_df = pd.DataFrame(columns=['gene', 'pos', 'Sanger'])
for k, v in sanger_muts_protease.items():
    dh = {'gene': 'protease', 'pos': k, 'Sanger': v}
    sanger_df = sanger_df.append(dh, ignore_index=True)
for k, v in sanger_muts_RT.items():
    dh = {'gene': 'RT', 'pos': k, 'Sanger': v}
    sanger_df = sanger_df.append(dh, ignore_index=True)

# MinVar results
minvar_file = 'minvar_results_%s/annotated_mutations.csv' % sys.argv[1]
minvar_muts_protease = {}
minvar_muts_RT = {}
with open(minvar_file) as csvfile:
    mut_list = csv.DictReader(csvfile, delimiter=',')
    for r in mut_list:
        if r['gene'] not in ['protease', 'RT']:
            continue
        int_freq = int(round(float(r['freq']) * 100, 0))
        pos = int(r['pos'])
        if r['gene'] == 'protease':
            try:
                minvar_muts_protease[pos].append('%s_%d' % (r['mut'], int_freq))
            except KeyError:
                minvar_muts_protease[pos] = ['%s_%d' % (r['mut'], int_freq)]
        if r['gene'] == 'RT':
            try:
                minvar_muts_RT[pos].append('%s_%d' % (r['mut'], int_freq))
            except KeyError:
                minvar_muts_RT[pos] = ['%s_%d' % (r['mut'], int_freq)]
minvar_df = pd.DataFrame(columns=['gene', 'pos', 'MinVar'])
for k, v in minvar_muts_protease.items():
    dh = {'gene': 'protease', 'pos': k, 'MinVar': v}
    minvar_df = minvar_df.append(dh, ignore_index=True)
for k, v in minvar_muts_RT.items():
    dh = {'gene': 'RT', 'pos': k, 'MinVar': v}
    minvar_df = minvar_df.append(dh, ignore_index=True)

# consensus B wild type sequences
wt = pd.DataFrame(columns=['pos', 'gene', 'wt'])
s = list(SeqIO.parse('protease.faa', 'fasta'))[0]
for i, a in enumerate(s):
    dh = {'pos': i + 1, 'gene': 'protease', 'wt': a}
    wt = wt.append(dh, ignore_index=True)

s = list(SeqIO.parse('RT.faa', 'fasta'))[0]
for i, a in enumerate(s):
    dh = {'pos': i + 1, 'gene': 'RT', 'wt': a}
    wt = wt.append(dh, ignore_index=True)

all_df = sanger_df.merge(minvar_df, how='outer', on=['gene', 'pos'])
all_df = all_df.merge(wt, how='inner', on=['gene', 'pos'])
all_df = all_df[all_df.pos <= 330]
all_df = all_df.sort(columns=['gene', 'pos'], ascending=[0, 1])
all_df['pos'] = all_df['pos'].astype(int)

all_df.to_csv(sys.stdout, sep='\t', index=False)
