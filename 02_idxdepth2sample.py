#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   02_idxdepth2sample.py
@Time    :   2021/04/28 19:32:37
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   inpux idxdepth putput dic; it's script for make a sample.txt for paragraph  input -i 
for example python 02_idxdepth2sample.py input_list samples.txt
'''

import datetime, sys, ast ,re
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

bamdir='03_bam/'

infile=sys.argv[1]
outfile=sys.argv[2]
f_in=open(infile,'r')
f_out=open(outfile,'w')

header=['id','path','depth','read_length']
f_out.write('\t'.join(header)+'\n')


for  i in f_in:
    f_tem = open(i.strip(),'r')
    id=re.split('[/]',i.strip())[1]
    tem=ast.literal_eval(f_tem.readlines()[0])['contigs']
    num=0
    for i in tem:
        num += float(i['depth'])
    num=num/12
    f_out.write(id+'\t'+bamdir+id+'.bam'+'\t'+str(round(num))+'\t'+'150'+'\n')
    print('Finished sample name : '+'\t'+id)

f_in.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')