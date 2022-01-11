export PATH=/public/agis/huangsanwen_group/chenglin/softwares/bwa-0.7.17:$PATH
export PATH=/vol1/agis/huangsanwen_group/lihongbo/software/bcftools-1.10.2/:$PATH
ref=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/DM_v6.1_all_chr.fa
#snp_vcf=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/44.snp.filt.recode.vcf.gz
name=DM

for i in {02..12};do
#gunzip 44pt.chr$i.pass.snp.vcf.gz
bcftools index 44pt.chr$i.pass.snp.vcf.gz
done
chenglin
##merge vcf
gatk GatherVcfs -I 44pt.chr01.1.pass.snp.vcf.gz -I 44pt.chr01.2.pass.snp.vcf.gz -I 44pt.chr02.pass.snp.vcf.gz -I 44pt.chr03.pass.snp.vcf.gz -I 44pt.chr04.pass.snp.vcf.gz -I 44pt.chr05.pass.snp.vcf.gz -I 44pt.chr06.pass.snp.vcf.gz -I 44pt.chr07.pass.snp.vcf.gz -I 44pt.chr08.pass.snp.vcf.gz -I 44pt.chr09.pass.snp.vcf.gz -I 44pt.chr10.pass.snp.vcf.gz -I 44pt.chr11.pass.snp.vcf.gz -I 44pt.chr12.pass.snp.vcf.gz -O 44.snp.vcf.gz
##filtering vcf soyebean 
vcftools --gzvcf $snp_vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --minDP 3 --recode-INFO-all --out 44.snp.filt
#--max-missing 过滤掉缺失率大于50%的位点
#--minQ 过滤掉低于30的质量的reads；
#--Mac 次要等位基因深度为3，过滤小于3的位点；
#--minDP 最低的深度

##make graph liner
ref=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/DM_v6.1_all_chr.fa
name=all44.liner
vg construct -r $ref -t 100  > $name.vg
##convert vg to xg
vg index -x $name.xg -g $name.gcsa -k 16 -t 100 -b /vol1/agis/huangsanwen_group/chenglin $name.vg
## for primary (liner reference) DM.paths is primary liner reference
vg index -T -G $name.gbwt $name.xg 
#vg gbwt -g $name.gg -x $name.xg $name.gbwt
vg gbwt -p -g $name.gg -o $name.gbwt -E -x $name.xg
#vg minimizer -t 16 -p -i $name.min -g $name.gbwt -G $name.gg
vg snarls $name.xg -T > $name.snarls
vg index -j $name.dist -s $name.snarls $name.xg 

##make graph second for 01_ref
vg construct -r $ref -t 20  > x.$name.vg
##index  -P option  will remove graph name
vg gbwt -p -g output.$name.gg -o output.$name.gbwt -E -x x.DM.xg
# vg gbwt -p -g output.gg -o output.gbwt -E -x x.xg -P 
vg minimizer -t 16 -p -i $name.min -g $name.gbwt -G $name.gg

#tell if there are sample names and #list the names
vg gbwt -M $name.gbwt
vg gbwt -S -L $name.gbwt

##simulate reads
name=all44.liner
f1=02_data/PG5002/s5002_HGVCKALXX_L1_1.clean.fq.gz
f2=02_data/PG5002/s5002_HGVCKALXX_L1_2.clean.fq.gz
vg sim -r -n 1000000 -a -s 12345 -p 570 -v 165 -i 0.00029 -x $name.xg -g $name.gbwt --sample-name ref \
--ploidy-regex "JTFH.*:0,KN70.*:0,Y:0,chrY_.*:0,chrEBV:0,.*:2" -F $f1 -F $f2 | vg annotate -p -x $name.xg -a - > sim.liner.gam

echo convert to fastq
vg view -X -a sim.liner.gam | gzip > sim.liner.fq.gz

echo format true position information
vg view -a sim.liner.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.liner.pos

#map
vg giraffe -x $name.xg -H $name.gbwt -m $name.min -G sim.liner.gam -i -t 22 >mapped.$name.gam 
vg gamcompare -r 100 -s <(vg annotate -m -x $name.xg -a mapped.$name.gam) sim.liner.gam 2>count | vg view -aj - > compared.$name.json
CORRECT_COUNT="$(sed -n '1p' count | sed 's/[^0-9]//g')"
SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"
MAPQ="$(grep mapping_quality\":\ 60 compared.$name.json | wc -l)"
MAPQ60="$(grep -v correctly_mapped compared.$name.json | grep mapping_quality\":\ 60 | wc -l)"
IDENTITY="$(jq '.identity' compared.$name.json | awk '{sum+=$1} END {print sum/NR}')"
GRAPH=potato
GBWT=liner-ref
READS=illumina
PARAM_PRESET=default
PAIRING=paired
SPEED=null
printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tscore\n" > report.$name.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats.$name.tsv
echo ${GRAPH} ${GBWT} ${READS} ${PAIRING} ${SPEED} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE}
printf "${GRAPH}\t${GBWT}\t${READS}\t${PAIRING}\t${SPEED}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${SCORE}\n" >> report.$name.tsv
jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.$name.json | sed 's/,/\t/g' | sed "s/$/\tgiraffe_primary_${GRAPH}${GBWT}${READS}${PAIRING}/" >> roc_stats.$name.tsv
grep -v 'null' roc_stats.$name.tsv > roc_stats.$name.ft.tsv



##bwa mapping 
ref=DM_v6.1_all_chr.fa
THREADS=20
READS=illumina
PAIRING=paired
SPEED=null
GRAPH=DM
ALGORITHM=bwa
GRAPH_NAME=DM_v6.1
#GRAPH_NAME=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/00_test/01_paragraph/DM_v6.1_all_chr.fa
bwa mem -t ${THREADS} -p ${ref} sim.liner.fq.gz > mapped.DM.raw.bam
samtools view -F 2048 -b mapped.DM.raw.bam > mapped.DM.bam
#bam=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/00_test/01_paragraph/03_bam/PG5002.bam 
vg inject -x $name.xg mapped.DM.bam -t 20 > mapped.DM.gam
vg view -aj mapped.DM.gam | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -aGJ - | vg annotate -m -x $name.xg -a - | vg gamcompare -r 100 -s - sim.liner.gam 2> count | vg view -aj - > compared.bwa.json

CORRECT_COUNT="$(grep correctly_mapped compared.bwa.json | wc -l)"
SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"
MAPQ="$(grep mapping_quality\":\ 60 compared.bwa.json | wc -l)"
MAPQ60="$(grep -v correctly_mapped compared.bwa.json | grep mapping_quality\":\ 60 | wc -l)"
IDENTITY="$(jq '.identity' compared.bwa.json | awk '{sum+=$1} END {print sum/NR}')"
echo ${GRAPH} ${READS} ${PAIRING} ${SPEED} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${SCORE}
printf "${GRAPH}\t${ALGORITHM}\t${READS}\t${PAIRING}\t-\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report_${ALGORITHM}.tsv
jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.bwa.json | sed 's/,/\t/g' | sed "s/$/\t${ALGORITHM}_${GRAPH}${READS}${PAIRING}/" | sed 's/single//g ; s/paired/-pe/g ; s/null/0/g' >> roc_stats_${ALGORITHM}.tsv
grep -v 'null' roc_stats_${ALGORITHM}.tsv > roc_stats_${ALGORITHM}.ft.tsv


##merge for roc figure
cat roc_stats_giraffe.tsv roc_stats.DM.paths.ft.tsv roc_stats_bwa.ft.tsv >roc-all.tsv


Rscript=/public/agis/huangsanwen_group/chenglin/softwares/miniconda3/envs/r42/bin/Rscript
$Rscript plot-qq.R roc_stats_bwa.ft.tsv roc_bwa.ft.qq.pdf
$Rscript plot-roc.R roc_stats_bwa.ft.tsv roc_bwa.ft.roc.pdf


