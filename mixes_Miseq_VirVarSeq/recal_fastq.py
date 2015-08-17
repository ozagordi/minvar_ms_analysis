#!/usr/local/bin/python3.4
import sys
import os
import warnings
import subprocess

def recalibrate_qualities(ref_file, bamfile, platform="ILLUMINA"):
    '''Invoke GATK BaseRecalibrator, also calling some high confidence variants.
    Follows vipr by Andreas Wilm'''

    bamstem, bamext = os.path.splitext(bamfile)
    assert bamext == '.bam', bamext
    recalfile = '%s_recal.bam' % bamstem
    if os.path.exists(recalfile):
        warnings.warn('file %s exists, not overwriting it' % recalfile)
        return recalfile
    # count mapped reads
    cml = "samtools view -c -F 4 %s" % bamfile
    mapped_reads = int(subprocess.check_output(cml, shell=True))

    # sample 10k reads and run lofreq to detect high confidence variants
    fraction_needed = min(1.0, round(float(10000) / mapped_reads, 3))
    cml = 'samtools view -s %f -F 4 -h -b -o subsampled.bam %s' % (fraction_needed, bamfile)
    subprocess.call(cml, shell=True)
    cml = 'samtools index subsampled.bam'
    subprocess.call(cml, shell=True)

    cml = 'lofreq call-parallel --pp-threads 6 -f %s -o tmp.vcf subsampled.bam' % ref_file
    subprocess.call(cml, shell=True)
    cml = 'lofreq filter -v 200 -V 2000 -a 0.40 -i tmp.vcf -o known.vcf'
    subprocess.call(cml, shell=True)
    os.remove('tmp.vcf')

    # need to add group information to reads
    print("@RG\tID:minvar\tSM:haploid\tLB:ga\tPL:%s" % platform, file=open('rg.txt', 'w'))
    cml = "samtools view -h %s |  cat rg.txt - | " % bamfile
    cml += "awk '{ if (substr($1,1,1)==\"@\") print; else printf \"%s\\tRG:Z:minvar\\n\",$0; }' | "
    cml += "samtools view -u - > grouped.bam"
    subprocess.call(cml, shell=True)
    os.remove('rg.txt')
    cml = 'samtools index grouped.bam'
    subprocess.call(cml, shell=True)

    refstem = os.path.splitext(ref_file)[0]
    cml = 'java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary R=%s O=%s.dict' % (ref_file, refstem)
    subprocess.call(cml, shell=True)

    # first pass to parse features
    cml = 'java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator  --maximum_cycle_value 600'
    cml += ' -I grouped.bam -l INFO -R %s -o recal_data.grp -knownSites known.vcf' % ref_file
    subprocess.call(cml, shell=True)

    # second pass to recalibrate
    cml = 'java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T PrintReads'
    cml += ' -R %s -I grouped.bam -BQSR recal_data.grp -o %s' % (ref_file, recalfile)
    subprocess.call(cml, shell=True)

    return recalfile


def main():
    recal_file = recalibrate_qualities('/home/ozagordi/Dropbox/Software/MinVar/minvar/db/consensus_B.fna', sys.argv[1])
    cml = 'samtools bam2fq -n %s' % recal_file
    subprocess.call(cml, shell=True)


if __name__ == '__main__':
    main()
