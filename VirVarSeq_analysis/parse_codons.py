#!/usr/local/bin/python3.4

import sys
import csv
import pandas as pd

def annogene(pdrow):
#    print(pdrow)
    pos = pdrow.POSITION
    if pos <= 99:
        gene = 'protease'
    elif 99 < pos and pos <= 539:
        gene = 'RT'
        pos -= 99
    elif 539 < pos and pos <= 659:
        gene = 'RNase'
        pos -= 539
    else:
        gene = 'Integrase'
        pos -= 659
    return pd.Series(dict(gene=gene, pos=pos))


infile, outfile = sys.argv[1:3]

codons = pd.read_csv(infile, sep='\t', header=0)
# discard codons at freq < 0.5% as recommended and sum synonymous mutations
codons = codons[codons.FREQ >= 0.005]
cfreq = codons.groupby(['POSITION', 'REF_AA', 'AA']).sum()
cfreq.reset_index(inplace=True)
# discard conserved mutations
cfreq = cfreq[cfreq.REF_AA != cfreq.AA]
# keep only relevant columns, sort and write
cfreq = cfreq.iloc[:, [0, 1, 2, -1]]
cfreq = cfreq.sort(columns=['POSITION', 'FREQ'], ascending=[1, 0])

# annotate adding gene/position info, based on region starting at protease 1
# 1 -> protease 1
# 99 -> protease 99
# 100 -> RT 1
# 539 -> RT 440
# 540 -> RNase 1
# 659 -> Rnase 120
# 660 -> Integrase 1
# 948 -> Integrase 289

annofreq = pd.concat([cfreq, cfreq.apply(annogene, axis=1)], axis=1)
#print(annofreq)
annofreq.to_csv(outfile, sep=',', float_format='%6.4f', index=False)
