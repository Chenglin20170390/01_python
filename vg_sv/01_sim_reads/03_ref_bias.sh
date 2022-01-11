#for test chr12; make a snp indel index for chr12
#workdir    /vol3/agis/huangsanwen_group/chenglin/work/6_pangenome/06_td_snp
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/gatk-4.1.9.0:$PATH
export PATH=/vol1/agis/huangsanwen_group/lihongbo/software/bcftools-1.10.2/:$PATH
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/bwa-0.7.17:$PATH

ref=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/DM_v6.1_all_chr.fa
vcf=44.snp.indel.sorted.vcf.gz
name=44.all

name=chr12-2
vcf=chr12.sort.filt2.vcf.gz
ref=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/chr12.fa

#combine vcf from chr / snp /indel
gunzip name
for i in 01.2 {02..11};do
bcftools index 44pt.chr${i}.pass.indel.vcf.gz
done
gatk GatherVcfs -I 44pt.chr01.1.pass.indel.vcf.gz -I 44pt.chr01.2.pass.indel.vcf.gz -I  44pt.chr02.pass.indel.vcf.gz \
-I 44pt.chr03.pass.indel.vcf.gz -I 44pt.chr04.pass.indel.vcf.gz -I 44pt.chr05.pass.indel.vcf.gz \
-I 44pt.chr06.pass.indel.vcf.gz -I 44pt.chr07.pass.indel.vcf.gz -I 44pt.chr08.pass.indel.vcf.gz \
-I 44pt.chr09.pass.indel.vcf.gz -I 44pt.chr10.pass.indel.vcf.gz -I 44pt.chr11.pass.indel.vcf.gz \
-I 44pt.chr12.pass.indel.vcf.gz -O 44.indel.vcf.gz

gatk GatherVcfs -I 44.snp.vcf.gz -I 44.indel.vcf.gz -O 44.snp.indel.vcf.gz
#gatk GatherVcfs -I 44pt.chr12.pass.snp.vcf.gz -I 44pt.chr12.pass.indel.vcf.gz -O chr12.all.vcf.gz
gatk SortVcf -I name -O name.out

##filtering vcf soyebean 
vcftools --gzvcf 44.snp.indel.vcf.gz --minQ 30 --recode --minDP 3 --recode-INFO-all --out 44.snp.indel
bgzip 44.snp.indel.vcf.gz

bcftools tabix $vcf


#vg construct
vg construct -r $ref -t 100 -v $vcf -a > x.$name.vg
#vg gbwt -p -g output.gg -o output.noP.gbwt -E -x x.xg
vg index -x x.$name.xg -k 16 -t 100 -b $tem x.$name.vg

##gbwt
##anthor recommend -n no -E; -n is meaning artifical reduce haplotype to 64 to downsample
vg gbwt -p -g output.$name.gg -o output.$name.gbwt -v $vcf -x x.$name.xg
#vg gbwt -p -g output.$name.gg -o output.$name.gbwt -x x.$name.xg -n 64 -l output.$name.gbwt
vg gbwt -p -g output.$name.gg -o output.$name.gbwt -x x.$name.xg -n 64 -l output.$name.gbwt
##nthor recommend
#vg gbwt -p -g output.$name.gg -o output.$name.gbwt -E -x x.$name.xg 
vg minimizer -t 20 -p -i output.$name.min -g output.$name.gbwt -G output.$name.gg
vg snarls x.$name.xg >x.$name.snarls
vg index -j $name.dist -s x.$name.snarls x.$name.vg 

#tell if there are sample names and #list the names

f1=02_data/PG5002/s5002_HGVCKALXX_L1_1.clean.fq.gz
f2=02_data/PG5002/s5002_HGVCKALXX_L1_2.clean.fq.gz

#giraffe
vg giraffe -H output.$name.gbwt -g output.$name.gg -m output.$name.min -d $name.dist -f $f1 -f $f2 -t 24 | vg surject -x x.$name.xg -b -t 24 - > mapped.giraffe.bam
samtools sort mapped.giraffe.bam -@ 20 > mapped.sorted.giraffe.bam
samtools index -@ 20 mapped.sorted.giraffe.bam
#bwa
bwa mem -t 24 -p $ref $f1 $f2 > mapped.bwa.bam
samtools sort mapped.bwa.bam -@ 20 > mapped.sorted.bwa.bam
samtools index -@ 20 mapped.sorted.bwa.bam

bcftools mpileup --threads 20 -f $ref -E -a DP -a SP -a ADF -a ADR -a AD  -O u  mapped.sorted.giraffe.bam | bcftools call --threads 20 -m -o calls.giraffe.vcf.gz -O z -v
bcftools mpileup --threads 20 -f $ref -E -a DP -a SP -a ADF -a ADR -a AD  -O u  mapped.sorted.bwa.bam | bcftools call --threads 20 -m -o calls.bwa.vcf.gz -O z -v

bcftools sort -O z calls.bwa.vcf.gz >calls.bwa.sorted.vcf.gz
bcftools sort -O z calls.giraffe.vcf.gz >calls.gif.sorted.vcf.gz

bcftools index -f calls.bwa.sorted.vcf.gz
bcftools index -f calls.gif.sorted.vcf.gz

# vi  rename bam name to mapped.sorted.bam
bcftools merge calls.gif.sorted.vcf.gz calls.bwa.sorted.vcf.gz -O z > all_calls.vcf.gz

##需要解压vcf 作图
gunzip all_calls.vcf.gz
python3 plot_allele_bia.py  all_calls.vcf allel_bias.svg


##toil vg
#source /public/home/chenglin/toilvenv/bin/activate
source /public/home/chenglin/toilvenv-36/bin/activate

toil-vg construct \
  --fasta chr12.fa   \
  --vcf chr12.sort.filt2.vcf.gz \
  --vcf_phasing chr12.sort.filt2.vcf.gz \
  --fasta_regions \
  --alt_paths \
  --out_name chr12 \
  --xg_index --gbwt_index --gcsa_index --trivial_snarls_index --distance_index \
  --gbwt_prune \
  --force_phasing True \
  --pangenome \
  --merge_graphs \
  --keep_vcfs \
  vg-test ./

##nthor recommend
vg gbwt -p -g output.$name.gg -o output.$name.gbwt -v $vcf -x x.$name.xg
vg gbwt -p -g output.$name.gg  -x x.$name.xg  output.$name.gbwt

#vg gbwt -p -g output.$name.gg -o output.$name.gbwt -E -x x.$name.xg 
vg minimizer -t 20 -p -i output.$name.min -g output.$name.gbwt -G output.$name.gg
vg snarls x.$name.xg >x.$name.snarls
vg index -j $name.dist -s x.$name.snarls x.$name.vg 



##new vg version
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/vg:$PATH
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/gatk-4.1.9.0:$PATH
export PATH=/vol1/agis/huangsanwen_group/lihongbo/software/bcftools-1.10.2/:$PATH
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/bwa-0.7.17:$PATH

name=chr12
ref=chr12.fa
vcf=chr12.sort.filt2.vcf.gz
tmp=/vol3/agis/huangsanwen_group/chenglin/work

vg construct -r $ref -t 20 -v $vcf -a > $name.vg
vg index -t 200 -p -L -b $tmp  -x $name.xg $name.vg 
#vg index -t 200 -p -b $tmp  -x $name.xg $name.vg 


vg snarls $name.xg >$name.snarls
vg index -p -L -b $tmp  $name.vg -s $name.snarls
vg gbwt -p -g $name.gg -o $name.gbwt -v $vcf -x $name.xg 

vg gbwt -p -g $name.gg -o $name.gbwt -x $name.xg -n 64 -l $name.gbwt
vg minimizer -t 20 -p -i $name.min -d $name.dist -g $name.gbwt -G $name.gg

f1=02_data/PG5002/s5002_HGVCKALXX_L1_1.clean.fq.gz
f2=02_data/PG5002/s5002_HGVCKALXX_L1_2.clean.fq.gz

#giraffe
vg giraffe -H $name.gbwt -g $name.gg -m $name.min -d $name.dist -f $f1 -f $f2 -t 24 | vg surject -x $name.xg -b -t 24 - > mapped.giraffe.bam
samtools sort mapped.giraffe.bam -@ 20 > mapped.sorted.giraffe.bam
samtools index -@ 20 mapped.sorted.giraffe.bam


##test

.gbwt
vg gbwt -p -g $name.gg -o $name-index.gbwt -x $name.xg -n 64 -l $name.gbwt
vg snarls $name.xg >$name.snarls
vg index -j $name.dist -s $name.snarls $name.vg
vg minimizer -t 20 -p -i $name.min -d $name.dist -g $name.gbwt -G $name.gg


##05_vg
vg gbwt -p -g $name.gg -o $name.gbwt -v $vcf -x $name.xg
vg snarls $name.xg >$name.snarls
vg index -j $name.dist -s $name.snarls $name.vg
vg minimizer -t 20 -p -i $name.min -d $name.dist -g $name.gbwt -G $name.gg
vg giraffe -H $name.gbwt -g $name.gg -m $name.min -d $name.dist -f $f1 -f $f2 -t 24 | vg surject -x $name.xg -b -t 24 - > mapped.giraffe.bam
warning[vg::giraffe]: Encountered 100000 ambiguously-paired reads before finding enough
                      unambiguously-paired reads to learn fragment length distribution. Are you sure
                      your reads are paired and your graph is not a hairball?
warning[vg::giraffe]: Finalizing fragment length distribution before reaching maximum sample size
                      mapped 661 reads single ended with 100000 pairs of reads left unmapped
                      mean: 0, stdev: 1
warning[vg::giraffe]: Cannot cluster reads with a fragment distance smaller than read distance
                      Fragment length distribution: mean=0, stdev=1
                      Fragment distance limit: 2, read distance limit: 200
warning[vg::giraffe]: Falling back on single-end mapping

samtools sort mapped.giraffe.bam -@ 20 > mapped.sorted.giraffe.bam
[E::hts_hopen] Failed to open file mapped.giraffe.bam
[E::hts_open_format] Failed to open file "mapped.giraffe.bam" : Exec format error

