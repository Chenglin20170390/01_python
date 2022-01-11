#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   01_clade_pick_from_Orthgroup.py
@Time    :   2021/06/15 19:27:41
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   提取特意clade,比如clade1.lst 里面C514, C356，C447，就是这三个材料的genecount其中任何一个不为0，而在除他们三以外的其它材料都为0的OG提出来
'''

import datetime, sys
from typing import ClassVar
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile=sys.argv[1]
outfile=sys.argv[2]
f_in=open(infile,'r')
f_out=open(outfile,'w')

orthgroup=open('01_clade_classification/Orthogroups.GeneCount.tsv','r')
list=[]
for i in f_in:
    list.append(i)

header=orthgroup.readline().strip().split()
f_out.write('\t'.join(header)+'\n')
list_idx=[]
for i in list:
    i=i.strip().split()[0]
    idx=header.index(i)
    list_idx.append(idx)
for line in orthgroup:
    line=line.strip().split()
    tol = int(line[-1])
    num=0
    sum=0
    for i in list_idx:
        if int(line[i]) != 0 :
            num += 1
            sum += int(line[i])
    if len(list_idx) == num and tol-sum == 0:
        f_out.write('\t'.join(line)+'\n')

print('Finished and output name :   ' + outfile)

f_in.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')



