#!/usr/local/bin/python3.4
'''Parses the output of Sierra (Stanford HIVdb algorithm) and puts it into
a more readable csv format
'''

import os
import csv
import sys
import shutil
import warnings
import re

filein = sys.argv[1]
md = {}
allpis = set([])
allrtis = set([])

with open(filein) as csvfile:
	mutreader = csv.DictReader(csvfile, delimiter='\t')
	for row in mutreader:
		sample = row['Seq ID']
		md[sample] = []
		rtis = row['NRTI'].split(',') + row['NNRTI'].split(',') + \
		row['RTI_OTHER'].split(',')
		rtis = [a for a in rtis if a != 'None']
		rtis = sorted(rtis, key=lambda mut: int(re.search('\w(\d*).*', mut).group(1)))

		pis = row['PI_MAJOR'].split(',') + row['PI_MINOR'].split(',') + \
		row['PI_OTHER'].split(',')
		pis = [a for a in pis if a != 'None']
		pis = sorted(pis, key=lambda mut: int(re.search('\w(\d*).*', mut).group(1)))

		for pi in pis:
			mtp = 'protease,%s' % pi
			#print(mtp)
			allpis.add(mtp)
			md[sample].append(mtp)
		for rti in rtis:
			mtp = 'RT,%s' % rti
			allrtis.add(mtp)
			#print(mtp)

			md[sample].append(mtp)

allpis = sorted(allpis, key=lambda mut: int(re.search('\w(\d*).*', mut.split(',')[1]).group(1)))
allrtis = sorted(allrtis, key=lambda mut: int(re.search('\w(\d*).*', mut.split(',')[1]).group(1)))
print('mix,gene,ref,pos,mut')
for k, v in md.items():
	for m in v:
		gene, mut = m.split(',')
		m_obj = re.search('(\w)(\d*)(.*)', mut)
		if m_obj.group(3).endswith('deletion'):
			print('%s,%s,%s,%s,%s' % (k, gene, m_obj.group(1), m_obj.group(2), '-'))
			continue
		else:
			for aa in m_obj.group(3):
				if m_obj.group(1) != aa:
					print('%s,%s,%s,%s,%s' % (k, gene, m_obj.group(1), m_obj.group(2), aa))

sys.exit()
for m in allpis:
	pres = [1 for x in md.keys() if m in md[x]]
	print(m + ',' + str(50*len(pres)))
	#print(md[sample]
for m in allrtis:
	pres = [1 for x in md.keys() if m in md[x]]
	print(m + ',' + str(50*len(pres)))
	#print(md[sample])
