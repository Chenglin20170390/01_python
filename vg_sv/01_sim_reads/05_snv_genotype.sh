export PATH=/public/agis/huangsanwen_group/chenglin/softwares/gatk-4.1.9.0:$PATH
export PATH=/vol1/agis/huangsanwen_group/lihongbo/software/bcftools-1.10.2/:$PATH
export PATH=/public/agis/huangsanwen_group/chenglin/softwares/bwa-0.7.17:$PATH


ref=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/DM_v6.1_all_chr.fa
vcf=44.snp.indel.ft.vcf.gz
name=all44
    Int NB_CRAM_CHUNKS = 8                              # Number of chunks to split the reads in
        Int MAX_CRAM_CHUNKS = 4                             # Number of chunks to actually analyze
        Boolean GIRAFFE_FAST_MODE = true                    # Set to 'false' to not use the fast mode for the read mapping with giraffe
        String GIRAFFE_RESCUE_MODE = "dozeu"                # Rescue mode for the giraffe mapper. Either "dozeu" (default) or "gssw" (slower but currently a bit more accurate)
        Boolean SNARL_CALLING = true                        # Return SV calls in the snarl context. Easier to merge results across samples, although requires 'bcftools norm -f REF.fa' to clean it up.
        Int CRAM_CONVERT_CORES = 4                          # Resources for the different tasks
        Int CRAM_CONVERT_DISK = 200
        Int MAP_CORES = 32
        Int MAP_DISK = 100
        Int MAP_MEM = 100
        Int MERGE_GAF_DISK = 400
        Int VGCALL_CORES = 16
        Int VGCALL_DISK = 200
        Int VGCALL_MEM = 100
        Int PREEMPTIBLE = 0                                # Number of attempts to use pre-emptible instances
        Boolean OUT_GAF = false    

f1=02_data/PG0027/s0027_HGLVTALXX_L4_1.clean.fq.gz
f2=02_data/PG0027/s0027_HGLVTALXX_L4_1.clean.fq.gz
sample_name=pg0027
thread=200
vcf=44.chr02.sorted.snp.indel.vcf.gz
name=chr02

##test for make snarls file
vg snarls -T $name.xg >$name.snarls
#generate gaf file for genotyping
vg giraffe -x $name.xg -m $name.min -d $name.dist -b fast --rescue-algorithm "dozeu" -N $sample_name --gbwt-name $name.giraffe.gbwt -C 500 -o gaf -f $f1 -f $f2 -t $thread | gzip > $sample_name.gaf.gz

# Create packed graph and genotype VCF
# GAF + GRAPH INDEXES + VCF_original -> PACK + VCF_genotyped
vg pack -x $name.xg -a $sample_name.gaf.gz -Q 5 -t $thread -o $sample_name.pack

vg call -k $sample_name.pack  -t $thread -s $sample_name --snarls $name.snarls -v $vcf $name.xg > $sample_name.gif.vcf
#vg call -k $sample_name.pack  -t $thread -s $sample_name  -v $vcf $name.xg > $sample_name.gif.vcf
bgzip $sample_name.gif.vcf
tabix -f -p vcf $sample_name.gif.vcv


##call at once (not genotype)
vg giraffe -x $name.xg -m $name.min -d $name.dist -b fast --rescue-algorithm "dozeu" -N $sample_name --gbwt-name $name.giraffe.gbwt -C 500 -o gam -f $f1 -f $f2 -t $thread > $sample_name.gam
vg convert $name.xg -p > $name.pg
# augment the graph, filtering out mappings and bases with quality < 5, and breakpoints with coverage < 4
vg augment $name.pg $sample_name.gam -m 2 -q 3 -Q 3 -t $thread -A $sample_name.aug.gam > $sample_name.aug.pg
# compute the snarls
vg snarls $sample_name.aug.pg > $sample_name.aug.snarls
# compute the support
vg pack -x $sample_name.aug.pg -g $sample_name.aug.gam -o $sample_name.aug.pack
# call variants on every path in the graph
vg call $sample_name.aug.pg -r $sample_name.aug.snarls -k $sample_name.aug.pack -s $sample_name > $sample_name.callsvcf 






#  vg mapping to sv graph to produce bam file
cd /vol3/agis/huangsanwen_group/chenglin/work/03_vg/00_ref_bias/03_vg_ref_bias
#for i in $(sed -n '13,20p' 402list);do
for i in $(sed -n '91,100p' 402list);do
name=all44
f1=`ls 02_data/$i/*1.clean.fq.gz` 
f2=`echo ${f1/1.clean/2.clean}`
outdir=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/00_ref_bias/04_gif_bam
sample_name=$i
nohup vg giraffe -H $name.giraffe.gbwt -g $name.gg -m $name.min -d $name.dist -x $name.xg -f $f1 -f $f2 -t 24 -o bam  | samtools sort -@ 24 > $outdir/$sample_name.sorted.gif.bam &
sleep 2s
#samtools index -@ 24 $outdir/$sample_name.sorted.gif.bam
done



##deepvariant
##variation calling snp and small indels
sample_name=PG0011
THREAD=100
mkdir -p /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/$sample_name
mkdir -p /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/tmp_gif/$sample_name
cd /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/tmp_gif/$sample_name
/public/software/singularity-3.5.2/bin/singularity exec \
-B /vol3/agis/huangsanwen_group/chenglin/work/03_vg/00_ref_bias/04_gif_bam/:/input  \
-B /vol3/agis/huangsanwen_group/chenglin/work/1_reference/:/ref \
-B /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/:/output \
  /public/software/singularity-3.5.2/docker-images/deepvariant-1.0.0.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/ref/DM_v6.1_all_chr.fa \
  --sample_name=$sample_name \
  --reads=/input/$sample_name.sorted.gif.bam \
  --output_vcf=/output/$sample_name/$sample_name.vcf.gz \
  --output_gvcf=/output/$sample_name/$sample_name.g.vcf.gz \
  --intermediate_results_dir /output/tmp_gif/$sample_name \
  --num_shards=${THREAD} \

#optional if you deepvariant not give option --sample_name=$sample_name you need reheader
#reheader
for i in PG0009;do
cd $i
bcftools view -h $i.vcf.gz >header.txt
sed -i 's/default/'$i'/g' header.txt 
bcftools reheader --threads 20 -h header.txt -o $i.rh.vcf.gz $i.vcf.gz
rm $i.vcf.gz
mv $i.rh.vcf.gz $i.vcf.gz
bcftools view -h $i.g.vcf.gz >header.txt
sed -i 's/default/'$i'/g' header.txt 
bcftools reheader --threads 20 -h header.txt -o $i.rh.g.vcf.gz $i.g.vcf.gz
rm $i.g.vcf.gz
mv $i.rh.g.vcf.gz $i.g.vcf.gz
rm header.txt
cd ..
done


##merge g.vcf to vcf by GLNeux
/public/software/singularity-3.5.2/bin/singularity exec \
-B /public/software/singularity-3.5.2/docker-images/glnexus_cli/:/opt/bin \
-B /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/merge/:/output \
/public/software/singularity-3.5.2/docker-images/deepvariant-1.0.0.sif /opt/bin/glnexus_cli \
--dir /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/tmp_merge \
--config DeepVariantWGS /vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf/*.g.vcf.gz | /vol1/agis/huangsanwen_group/lihongbo/software/bcftools-1.10.2/bcftools view - | bgzip -c  > 207.vcf.gz


## sample ti/tv ratio statistic
DIR=/vol3/agis/huangsanwen_group/chenglin/work/03_vg/05_snp/02_gif_vcf
for SAMPLE in $(cat list);do
  bcftools stats -f PASS \
    ${DIR}/${SAMPLE}.vcf.gz \
  > ${DIR}/${SAMPLE}.stats
done