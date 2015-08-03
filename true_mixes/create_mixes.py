#!/usr/local/bin/python3.4
import re
import sys
import csv

mixes = {}
mixes[9] = {'JRCSF': 0.52, 'YU2': 0.26, 'INP0223': 0.13, 'INP0224': 0.06, 'INP0157': 0.03}
mixes[10] = {'JRCSF': 0.20, 'YU2': 0.20, 'INP0223': 0.20, 'INP0224': 0.20, 'INP0157': 0.20}
mixes[11] = {'JRCSF': 0.60, 'YU2': 0.10, 'INP0223': 0.10, 'INP0224': 0.10, 'INP0157': 0.10}
mixes[12] = {'JRCSF': 0.80, 'YU2': 0.05, 'INP0223': 0.05, 'INP0224': 0.05, 'INP0157': 0.05}
mixes[13] = {'JRCSF': 0.90, 'YU2': 0.025, 'INP0223': 0.025, 'INP0224': 0.025, 'INP0157': 0.025}

mixes[17] = mixes[9]
mixes[18] = mixes[10]
mixes[19] = mixes[11]
mixes[20] = mixes[12]
mixes[21] = mixes[13]

mix = int(sys.argv[1])
freq = {}
mutfile = 'sequenced_clones.csv'
with open(mutfile) as csvfile:
        mutreader = csv.DictReader(csvfile, delimiter=',')
        for row in mutreader:
            clone = row['mix'].split('_')[-1]
            gene = row['gene']
            mut = '_'.join([row['gene'], row['ref'], str(row['pos']), row['mut']])
            freq[mut] = freq.get(mut, 0) + mixes[mix][clone]

print('gene,wt,pos,mut,freq')
for k, v in freq.items():
    # print(k, v)
    ksp = k.split('_')
    print(','.join(ksp) + ',%f' % round(v, 4))
