#!/usr/bin/env python3.4

import glob
import subprocess

dirs = glob.glob('minvar_results_*')
stats = {}
print('sample,subtype,percent_mapped')
for d in dirs:
    wc_out = subprocess.check_output('wc -l %s/high_quality.fastq' % d, shell=True)
    hq_reads = int(int(wc_out.split()[0])/4)
    idx_out = subprocess.check_output('samtools stats %s/hq_2_cons_sorted_recal.bam | grep ^SN | cut -f 2-' % d, shell=True).decode('utf-8')
    lines = idx_out.split('\n')
    for l in lines:
        try:
            stats[l.split(':')[0]] = l.split(':')[1].strip()
        except IndexError:
            pass

    sid = d.split('_')[-1]
    seqs = stats['sequences']
    mapped = stats['reads mapped']
    unmapped = stats['reads unmapped']
    #print('\t'.join([sid, pat_dict[sid], seqs, mapped, unmapped, str(100 * int(mapped) / int(seqs))]))
    print(','.join([sid, 'B', str(100 * int(mapped) / int(seqs))]))
