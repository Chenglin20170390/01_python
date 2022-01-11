#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   19_bin_fre_introgression.py
@Time    :   2021/06/03 19:28:02
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   None
'''

import datetime, sys
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile=sys.argv[1]
outfile=sys.argv[2]
f_in=open(infile,'r')
f_out=open(outfile,'w')

dic_po={'chr01':88591686,'chr02':46102915,"chr03":60707570,'chr04':69236331,'chr05':55599697,'chr06':59091578,'chr07':57639317,'chr08':59226000,'chr09':67600300,'chr10':61044151,'chr11':46777387,'chr12':59670755}
i=0
num=0
for line in f_in:
    if line.startswith('#'):
        continue
    if pos2 < 59670755:
        pos1=1000*i
        pos2=1000*(i+1)
        line=line.strip().split()[1]
        if pos1 < int(line) < pos2:
            num += 1
        elif 1000*(i+1) < int(line) < 1000*(i+2):
            f_out.write('chr12'+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(num)+'\n')
            i += 1
            num=1
        else:
            t=2
            for t < 1000:
                if 1000*(i+t) < int(line) < 1000*(i+t):
                    i=i+t
                    num=1
            
            
            



f_in.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')