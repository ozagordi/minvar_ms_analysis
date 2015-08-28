#!/usr/bin/env python3.4

from Bio import SeqIO


dna_code = {
    'A': set(['A']),
    'C': set(['C']),
    'G': set(['G']),
    'T': set(['T']),

    'R': set(['G', 'A']),
    'Y': set(['T', 'C']),
    'M': set(['A', 'C']),
    'K': set(['G', 'T']),
    'S': set(['G', 'C']),
    'W': set(['A', 'T']),

    'H': set(['A', 'C', 'T']),
    'B': set(['C', 'G', 'T']),
    'V': set(['A', 'C', 'G']),
    'D': set(['A', 'G', 'T']),
    'N': set(['A', 'C', 'G', 'T']),
    '-': set(['A', 'C', 'G', 'T'])
}
# 64 codons + '---'
translation_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*', '---': '-'}


def disentangle(codon):
    '''Codon with ambiguous bases gets expanded'''
    from itertools import product
    codons = []
    for x in product(''.join(dna_code[codon[0]]), ''.join(dna_code[codon[1]]), ''.join(dna_code[codon[2]])):
        codons.append(''.join(x))
    return(codons)

def parse_muts(s, gene):
    '''Input sequence with ambiguous nucleotides, output a list of
    mutations w.r.t. consensus B. If ambiguous nt includes the wild type,
    output both'''
    import os
    import subprocess
    from Bio import AlignIO

    cml = 'needle %s.fna asis:%s -gapopen 20.0 -gapextend 5.0 -outfile tmp.out -aformat3 fasta -auto' % (gene, s)
    subprocess.call(cml, shell=True)
    al = AlignIO.read('tmp.out', 'fasta')
    os.remove('tmp.out')
    all_muts = []
    if gene == 'protease':
        assert not al[0].seq.startswith('-')
        coords = range(0, 297, 3)
    elif gene == 'RT':
        coords = range(297, 1287, 3)
    for i in coords:
        wt = str(al[0][i:i + 3].seq.translate())
        cod = al[1][i:i + 3]
        pos = int(i / 3 + 1)
        if gene == 'RT':
            pos -= 99
        aa = str(cod.seq.translate())
        if pos == 67 and gene == 'RT':
            print(wt, aa, cod.seq)
        if aa == wt:
            continue
        # Biopython translates to ambiguous amminoacids, as in
        # RAT [AAT / GAT] -> B (aspartic acid or asparagine). We don't want it
        if aa == 'X' or aa not in translation_table.keys():
            am_cod = disentangle(str(cod.seq))
            to_add = set([','.join([gene, wt, str(pos), translation_table[c]]) for c in am_cod])
            freq = 1. / len(to_add)
            for t in to_add:
                all_muts.append(t + ',%4.3f' % freq)
                #print(','.join([gene, wt, str(pos), translation_table[c]]))
        else:
            to_add = ','.join([gene, wt, str(pos), aa, '1.0'])
            # if to_add not in all_muts:
            all_muts.append(to_add)
    return all_muts

mixes = ['Mix_%d' % i for i in [9, 11, 12, 13, 17, 18, 19, 20, 21]]
infile = 'Mix_1-24_CS_Aligned.fas'
seqs = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
for mix in mixes:
    s = seqs[mix]
    prot_muts = parse_muts(str(s.seq), 'protease')
    rt_muts = parse_muts(str(s.seq), 'RT')
    oh = open('%s_mutations.csv' % mix, 'w')
    oh.write('gene,wt,pos,mut,freq\n')
    for m in prot_muts:
        oh.write(m + '\n')
    for m in rt_muts:
        oh.write(m + '\n')
