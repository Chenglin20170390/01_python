#!/usr/bin/env python3

import sys,re

i1 = open(sys.argv[1])  #  80_a_stringtie.xls
i2 = open(sys.argv[2])  #  80_b_stringtie.xls
i3 = open(sys.argv[3])  #  80_c_stringtie.xls
i4 = open(sys.argv[4])  #  80_c_stringtie.xls
i5 = open(sys.argv[5])  #  80_c_stringtie.xls
i6 = open(sys.argv[6])  #  80_c_stringtie.xls
i7 = open(sys.argv[7])  #  80_c_stringtie.xls

d1 = {}
for line in i1:
	if '#' not in line:
		if line.split()[2] == "transcript":
			sline = line.strip().split()
			#re1 = re.search('gene_id \"(\w+)\"; ', line)
			re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
			#gene_id = re1.group(1)
			gene_id = sline[9][1:-2]
			fpkm = re2.group(1)
			d1[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d2 = {}
for line in i2:
	if '#' not in line and line.split()[2] == "transcript":
		sline = line.strip().split()
#		re1 = re.search('gene_id \"(\w+)\"; ', line)
		re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#		gene_id = re1.group(1)
		gene_id = sline[9][1:-2]
		fpkm = re2.group(1)
		d2[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d3 = {}
for line in i3:
	if '#' not in line and line.split()[2] == "transcript":
		sline = line.strip().split()
#		re1 = re.search('gene_id \"(\w+)\"; ', line)
		re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#		gene_id = re1.group(1)
		gene_id = sline[9][1:-2]
		fpkm = re2.group(1)
		d3[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d4 = {}
for line in i4:
    if '#' not in line and line.split()[2] == "transcript":
        sline = line.strip().split()
#       re1 = re.search('gene_id \"(\w+)\"; ', line)
        re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#       gene_id = re1.group(1)
        gene_id = sline[9][1:-2]
        fpkm = re2.group(1)
        d4[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d5 = {}
for line in i5:
    if '#' not in line and line.split()[2] == "transcript":
        sline = line.strip().split()
#       re1 = re.search('gene_id \"(\w+)\"; ', line)
        re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#       gene_id = re1.group(1)
        gene_id = sline[9][1:-2]
        fpkm = re2.group(1)
        d5[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d6 = {}
for line in i6:
    if '#' not in line and line.split()[2] == "transcript":
        sline = line.strip().split()
#       re1 = re.search('gene_id \"(\w+)\"; ', line)
        re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#       gene_id = re1.group(1)
        gene_id = sline[9][1:-2]
        fpkm = re2.group(1)
        d6[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

d7 = {}
for line in i7:
    if '#' not in line and line.split()[2] == "transcript":
        sline = line.strip().split()
#       re1 = re.search('gene_id \"(\w+)\"; ', line)
        re2 = re.search('FPKM \"(\d+\.\d+)\"; ', line)
#       gene_id = re1.group(1)
        gene_id = sline[9][1:-2]
        fpkm = re2.group(1)
        d7[gene_id+'~'+sline[0]+'~'+sline[3]+'~'+sline[4]] = float(fpkm)

print('#Gene\tChr\tStart\tEnd\tRoot\tStem\tLeaf\tMale.flower\tFemale.flower\tFruit\tTendril')
for k in d1:
	if k in d2 and k in d3 and k in d4 and k in d5 and k in d6 and k in d7:
		if 'pilon' in k:
			print(k.split('~')[0],k.split('~')[1],k.split('~')[2],k.split('~')[3],d1[k],d2[k],d3[k],'%.4f' % (float(d1[k]+d2[k]+d3[k])/3),sep='\t')
		else:
			print(k.split('~')[0],k.split('~')[1],k.split('~')[2],k.split('~')[3],d1[k],d2[k],d3[k],d4[k],d5[k],d6[k],d7[k],sep='\t')

