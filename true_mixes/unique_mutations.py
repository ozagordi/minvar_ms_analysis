#!/usr/bin/env python3.4

filein = 'sequenced_clones.csv'
counts = {}
clone = {}
for l in open(filein):
    if l.startswith('mix'):
        continue
    lsp = l.strip().split(',')
    clone_here, mut = lsp[0].split('_')[-1], ','.join(lsp[1:])
    counts[mut] = counts.get(mut, 0) + 1
    clone[mut] = clone_here

oh = open('unique_mutations.csv', 'w')
oh.write('clone,gene,ref,pos,mut\n')
for k, v in counts.items():
    if v == 1:
        print('%s,%s' % (clone[k], k), file=oh)
