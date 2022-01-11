#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
---------------------------------
   File Name    :   11_sniffles_vg.vcf
   Author       :   Administrator Lin Cheng
   Date         :   2020/12/16 0016
   Description  :   
---------------------------------
'''
import sys, os, re,datetime, glob

print("**********************************************************************")
print('Start Time: ' + str(datetime.datetime.now()))
print("**********************************************************************")

infile = sys.argv[1]
outfile = sys.argv[2]
f_in = open(infile, 'r')
f_out=open(outfile,'w')


def DNA_complement2(sequence):
    sequence=sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
ref='/vol3/agis/huangsanwen_group/chenglin/work/1_reference/DM_v6.1_all_chr.fa'
for line in f_in:
    if line.startswith('##'):
        f_out.write(line)
    elif line.startswith('#CHROM'):
         header_f=line.strip().split('\t')
         print(header_f)
         idx_all=header_f.index('FORMAT')
         for i in range (idx_all+1, len(header_f)):
             header.append(header_f[i])
         out_line='\t'.join(header)+'\n'
         #print(out_line)
         f_out.write(out_line)
    else:
        in_line=line.strip().split("\t")
        out_line=in_line[:9]   #DT
        for i in range(1, len(header_f)-idx_all):
            #print(len(header_f))
            #print(idx_all)
            #print(i)
            buf = in_line[idx_all+i].strip().split(":")[0]
            buf=str(buf)
            out_line.extend([buf])
        pos=re.split(';|=',in_line[7])
        idx_type=int(pos.index('SVTYPE'))+1
        sv_type=pos[idx_type]
        if sv_type!='TRA':
            idx=int(pos.index('END'))+1
            seq=''
            end=pos[idx]
            if sv_type=='DEL' and 0<=int(end)-int(out_line[1]):
                n=os.popen('samtools faidx '+ref+' '+out_line[0]+ ':'+out_line[1]+'-'+end).readlines()[1:]
                for i in n:
                    i=i.strip()
                    seq+=i
                seq=seq.strip()
                n1=seq[0] 
                out_line[4]=n1
                out_line[3]=seq
                out_line='\t'.join(out_line) + '\n' 
                f_out.write(out_line)
            if sv_type=='INS' and 0<=int(end)-int(out_line[1]):
                n=os.popen('samtools faidx '+ref+' '+out_line[0]+ ':'+out_line[1]+'-'+out_line[1]).readlines()[1].strip()
                if out_line[4]!='INS' and out_line[4]!='DEL' and out_line[4]!='INV' and out_line[4]!='DUP' and out_line[4]!='TRA':
                    out_line[3]=n
                    ##pick up first one #out_line[4]=out_line[4].strip().split(',')[0]
                    #out_line[4]=out_line[4]
                    out_line='\t'.join(out_line) + '\n'
                    f_out.write(out_line)
            if sv_type=='INV' and 0<= int(end) - int(out_line[1]):
                n=os.popen('samtools faidx '+ref+' '+out_line[0]+ ':'+out_line[1]+'-'+end).readlines()[1:]
                seq=''
                for i in n:
                    i=i.strip()
                    seq+=i
                seq=seq.strip()
                seq2=DNA_complement2(seq)
                out_line[3]=seq
                out_line[4]=seq2
                out_line='\t'.join(out_line) + '\n'
                f_out.write(out_line)
            if sv_type=='DUP' and 0<=int(end)-int(out_line[1]):
                n=os.popen('samtools faidx '+ref+' '+out_line[0]+ ':'+out_line[1]+'-'+end).readlines()[1:]
                seq=''
                for i in n:
                    i=i.strip()
                    seq+=i
                seq=seq.strip()
                out_line[3]=seq[0]
                out_line[4]=seq
                out_line='\t'.join(out_line) + '\n'
                f_out.write(out_line)

f_in.close()
f_out.close()
print("**********************************************************************")
print('End Time: ' + str(datetime.datetime.now()))
print("**********************************************************************")