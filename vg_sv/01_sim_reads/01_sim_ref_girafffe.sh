
#graph construct
vg construct -r DM_v6.1_all_chr.fa -v hsvlr.vcf.gz -t 20 > xcl.vg

##index  -P option  will remove graph name
vg gbwt -p -g output.gg -o output.noP.gbwt -E -x x.xg
# vg gbwt -p -g output.gg -o output.gbwt -E -x x.xg -P 
vg minimizer -t 16 -p -i output.min -g output.gbwt -G output.gg
##for alt index hgsvc
vg index  -t 20 -x output.alt.xg -L x.vg
vg prune x.vg --threads 32 --mapping output.pruned.mapping --unfold-paths --gbwt-name output.gbwt --progress > output.pruned.vg


#tell if there are sample names
vg gbwt -M x.gbwt
#list the names
vg gbwt -S -L x.gbwt

f1=02_data/PG5002/s5002_HGVCKALXX_L1_1.clean.fq.gz
f2=02_data/PG5002/s5002_HGVCKALXX_L1_2.clean.fq.gz
vg sim -r -n 10000 -a -s 12345 -p 570 -v 165 -i 0.00029 -x output.xg -g output.gbwt --sample-name ref \
--ploidy-regex "JTFH.*:0,KN70.*:0,Y:0,chrY_.*:0,chrEBV:0,.*:2" -F $f1 -F $f2 | vg annotate -p -x output.xg -a - > sim.gam

echo convert to fastq
vg view -X -a sim.gam | gzip > sim.fq.gz

echo format true position information
vg view -a sim.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.pos

## mapping 
vg giraffe -x output.xg -H output.gbwt -m output.min -G sim.gam  -b default -i > mapped.gam
vg giraffe -x output.xg -H output.gbwt -m output.min -G sim.gam  -b fast -i > mapped.gam
vg gamcompare -r 100 -s <(vg annotate -m -x output.xg -a mapped.gam) sim.gam 2>count | vg view -aj - > compared.json
CORRECT_COUNT="$(sed -n '1p' count | sed 's/[^0-9]//g')"
SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"
MAPQ="$(grep mapping_quality\":\ 60 compared.json | wc -l)"
MAPQ60="$(grep -v correctly_mapped compared.json | grep mapping_quality\":\ 60 | wc -l)"
IDENTITY="$(jq '.identity' compared.json | awk '{sum+=$1} END {print sum/NR}')"
GRAPH=potato
GBWT=full
READS=illumina
PARAM_PRESET=default
PAIRING=paired
SPEED=null
printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats.tsv
echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${SPEED} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE}
printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${SPEED}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report.tsv
jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.json | sed 's/,/\t/g' | sed "s/$/\tgiraffe_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> roc_stats.tsv



for i in {02..12};do
nohup vcftools --gzvcf ../44pt.chr$i.pass.snp.vcf.gz --SNPdensity 1000000 --out chr$i.density &
done

vcftools --vcf chr01.snp.vcf --SNPdensity 1000000 --out chr01.density
cat chr* |grep 'chr' >chr.all.snp.density
python add.colum.py