# this was run with sdev -m 48g on scg

# load modules
#-------------------------------
ml load plink2
ml load bcftools
ml load gatk

# set params
#-------------------------------
ukbb=/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome
elite=/oak/stanford/groups/euan/projects/elite/bcbio/projects/project_n267_gatk_joint/elite_n267_joint.vcf.gz

# transform ukbb to vcf
#-------------------------------
# plink2  --export vcf-4.2 vcf-dosage=DS --out ukbb_exome --pfile /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome
## bgzip
# bgzip ukbb_exome.vcf
# Or use our sh file
#
## add index
# bcftools index --csi ukbb_exome.vcf.gz

# liftover elite's data
#-------------------------------
# add tbi index file
gatk IndexFeatureFile -F $elite -O $elite.tbi

# get the chain files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
wget https://raw.githubusercontent.com/broadgsa/gatk/master/public/chainFiles/b37tohg19.chain
wget https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain
gunzip hg19ToHg38.over.chain.gz
mkdir chain_files
mv hg19ToHg38.over.chain chain_files
mv b37tohg19.chain chain_files
mv b37ToHg38.over.chain chain_files

# set the target references
b19Ref=/oak/stanford/groups/mrivas/public_data/genomes/hg19/hg19.fa
b38Ref=/oak/stanford/groups/mrivas/public_data/genomes/hg38/hg38.fa
# add genome dictionaries
mkdir hg38_ref
cp $b38Ref hg38_ref/
gatk CreateSequenceDictionary -R hg38_ref/hg38.fa > hg38_ref/create_dict.log
mkdir hg19_ref
cp $b19Ref hg19_ref/
gatk CreateSequenceDictionary -R hg19_ref/hg19.fa > hg19_ref/create_dict.log

# Run liftover1: grch37 to hg19
gatk --java-options "-Xmx24g" LiftoverVcf \
    -I $elite -O elite_n267_joint_hg19.vcf.gz -C chain_files/b37tohg19.chain \
    -R hg19_ref/hg19.fa --REJECT elite_grch37ToHg19_rejects.vcf.gz > liftover_grch37tohg19.log 2> liftover_grch37tohg19.err
# Run liftover2: hg19 to hg39
gatk --java-options "-Xmx24g" LiftoverVcf \
    -I elite_n267_joint_hg19.vcf.gz -O elite_n267_joint_hg38.vcf.gz \
    -C chain_files/hg19ToHg38.over.chain -R hg38_ref/hg38.fa \ 
    --REJECT elite_Hg19ToHg38_rejects.vcf.gz > liftover_hg19tohg38.log 2> liftover_hg19tohg38.err
# Directly from 37 to hg38
gatk --java-options "-Xmx40g" LiftoverVcf  -I $elite -O elite_n267_joint_hg38.v2.vcf.gz \
    -C chain_files/b37ToHg38.over.chain  -R hg38_ref/hg38.fa \ 
    --REJECT elite_grch37ToHg38_rejects.vcf.gz > liftover_grch37tohg38.log 2> liftover_grch37tohg38.err

# Comment: liftover can create variants with the same position and alleles, see example:
# (validated using the online tool)
# chr1:1582880-1582881
# chr1:1646113-1646114
# both are mapped to the same position:
# chr1:1714674-1714675

# preprocess elite's vcf
#-------------------------------

# add index
bcftools index --csi elite_n267_joint_hg38.vcf.gz

# filter the vcf using GATK's best practices + to have reasonable call rates
bcftools filter -e "QUAL<30 || ExcessHet>60 || MAC[0]<2 || AN < 100 || QD<2 || FS>60 || SOR>3 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8" elite_n267_joint_hg38.vcf.gz -O z -o elite_n267_joint_hg38.filtered.vcf.gz
bcftools index --csi elite_n267_joint_hg38.filtered.vcf.gz

## preprocess ukbb's vcf
##-------------------------------
## filter the vcf
#bcftools filter -i "MAC[0]>4" -i "AN>10000" ukbb_exome.vcf.gz -O z -o ukbb_exome.filtered.vcf.gz
#bcftools index --csi ukbb_exome.filtered.vcf.gz

# elite - remove multiallelic variants
bcftools view -Oz --max-alleles 2 elite_n267_joint_hg38.filtered.vcf.gz > elite_n267_joint_hg38.filtered.biallelic.vcf.gz
bcftools index --csi elite_n267_joint_hg38.filtered.biallelic.vcf.gz
# rm elite_n267_joint_hg38.filtered.biallelic.vcf.gz.csi

# elite - add variant ids
bcftools annotate -O z --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' elite_n267_joint_hg38.filtered.biallelic.vcf.gz > elite_n267_joint_hg38.filtered.biallelic.anno.vcf.gz
bcftools index --csi elite_n267_joint_hg38.filtered.biallelic.anno.vcf.gz

# Analyze the vcfs
#-------------------------------
# F1=elite_n267_joint_hg38.filtered.biallelic.anno.vcf.gz
# F2=ukbb_exome.vcf.gz

# Slow and has an empty output:
# bcftools isec -n +2 $F1 $F2 -p . -O z
# bcftools merge 0000.vcf.gz 0001.vcf.gz -O z -o merged.vcf.gz
# plink2 --vcf merged.vcf.gz --no-sex --maf 0.01 --make-bed --out elite.ukbb.merged
# rm 0001.vcf.gz
# rm 0000.vcf.gz
# plink2 --bfile elite.ukbb.merged --indep-pairwise 500 10 0.3 --maf 0.01 --out elite.ukbb.merged
# plink2 --bfile elite.ukbb.merged --extract elite.ukbb.merged.pruned.prune.in --make-bed --pca --out elite.ukbb.merged.pruned



