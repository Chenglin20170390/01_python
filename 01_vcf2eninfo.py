#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   01_vcf2eninfo.py
@Time    :   2021/04/28 09:55:25
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   for convert sv_vcf (ref/alt) lenght to add to info file if the info not including(END=？)
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


for line in f_in:
    if line.startswith("#"):
        f_out.write(line)
    else:
        line=line.strip().split()
        ref=line[3]
        alt=line[4]
        if len(ref) == 1 and len(alt) > 1:
            svtype='SVTYPE=INS'
            le_n = 'SVLEN='+str(len(alt))
            end = 'END='+str(int(line[1])+len(alt)-1)
            line[7]=svtype+';'+le_n+';'+end
            f_out.write('\t'.join(line)+'\n')
        if len(ref) > 1 and len(alt) == 1:
            svtype='SVTYPE=DEL'
            le_n = 'SVLEN='+str(len(ref))
            end = 'END='+str(int(line[1])+len(ref)-1)
            line[7]=svtype+';'+le_n+';'+end
            f_out.write('\t'.join(line)+'\n')
        if len(ref) == len(alt) :
            svtype='SVTYPE=INV'
            le_n = 'SVLEN='+str(len(alt))
            end = 'END='+str(int(line[1])+len(alt)-1)
            line[7]=svtype+';'+le_n+';'+end
            f_out.write('\t'.join(line)+'\n')
        if len(ref) != 1 and len(alt) != 1 and  len(alt) % len(ref) ==0:
            svtype='SVTYPE=DUP'
            le_n = 'SVLEN='+str(len(alt))
            end = 'END='+str(int(line[1])+len(alt)-1)
            line[7]=svtype+';'+le_n+';'+end
            f_out.write('\t'.join(line)+'\n')
        
        

f_in.close()
f_out.close()

print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')