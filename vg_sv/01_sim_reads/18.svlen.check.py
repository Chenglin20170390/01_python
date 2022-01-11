#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   18.svlen.check.py
@Time    :   2021/05/31 11:17:59
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin20170390@163.com
@License :   (C)Copyright 2019-2020, CAAS ShenZhen
@Desc    :   for checking and filtering sv len from population file
'''

import datetime, sys,re 
# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile=sys.argv[1]
outfile=sys.argv[2]
f_in=open(infile,'r')
f_out=open(outfile,'w')

for line in f_in:
    if  line.startswith("#"): 
        f_out.write(line)
    else:
        line=line.strip().split()
        st = line[2]
        info = re.split(";=",line[9])
        ed_idx=info.index("END")
        ed=info[int(ed_idx)+1]
        if int(ed) >=int(st):
            f_out.write(line)
f_in.close()
f_out.close()

        
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')