

#  vg mapping to sv graph
f1=02_data/PG5002/s5002_HGVCKALXX_L1_1.clean.fq.gz
f2=02_data/PG5002/s5002_HGVCKALXX_L1_2.clean.fq.gz
vg giraffe -x output.xg -H output.gbwt -m output.min -t 20 -b default -f $f1 -f $f2 | vg surject -x output.xg -b -t 20 - > mapped.5002.bam
samtools sort mapped.5002.bam -@ 20 > mapped.5002.sorted.bam
samtools index -@ 20 mapped.5002.sorted.bam


##check coverage
samtools flagstat PG5002.bam > PG5002.bwa.sum
samtools flagstat mapped.5002.sorted.bam> PG5002.gif.sum

