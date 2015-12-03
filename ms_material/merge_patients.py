#!/usr/bin/env python3.4
import pandas as pd

for pat in range(1, 11):
    in_454 = '../patients_454/summary/comparison_%d.tsv' % pat
    in_miseq = '../patients_Miseq/summary/comparison_%d.tsv' % pat
    df_454 = pd.read_csv(in_454, sep='\t')
    df_miseq = pd.read_csv(in_miseq, sep='\t')

    full = pd.merge(df_454, df_miseq, how='outer', on=['gene', 'pos', 'wt'],
                    suffixes=('_454', '_Miseq'))
    full = full.drop('Sanger_454', axis=1)
    full = full.rename(columns={'Sanger_Miseq': 'Sanger'})
    full['pos'] = full['pos'].astype(int)
    full = full.sort(columns=['gene', 'pos'], ascending=[0, 1])
    # full = full.where(pd.notnull(full), None)
    full.to_csv('patient_%d_merged.tsv' % pat, sep='\t', index=False)
