# set params
#-------------------------------
elite=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_n267_joint_hg38.filtered.biallelic.anno.vcf.gz
out=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_analyze_vcf

# extract ref file
# slow - run once or use a prev calculated one
# less $elite | grep -v "#" | cut -f 3,4 > elite.ref.txt
# creating plink files rquires two steps
# 1. identify duplicated snp ids as a result from the liftover
# this command fails but create a file with the ids to be excluded
plink2 --vcf $elite  --make-bed --out ${out}/elite_exome --allow-extra-chr --rm-dup
# 2. Run the plink conversion with the reference and the exclusion of duplicated variant ids
plink2 --vcf $elite  --make-bed --out ${out}/elite_exome --allow-extra-chr --rm-dup --exclude elite_exome.rmdup.mismatch --ref-allele elite.ref.txt 
plink2 --vcf $elite  --make-pgen --out ${out}/elite_exome --allow-extra-chr --rm-dup --exclude elite_exome.rmdup.mismatch --ref-allele elite.ref.txt

# get relatedness estimates (May 2020: no pairs identified)
# first prune and maf filters
plink --bfile ${out}/elite_exome --maf 0.01 --indep-pairwise 500 10 0.1  --out ${out}/elite_exome --allow-extra-chr
plink --bfile ${out}/elite_exome --extract ${out}/elite_exome.prune.in --chr 1-22  --genome --min 0.2 --out ${out}/elite_exome --allow-extra-chr

# inbreeding coeffs
plink --bfile ${out}/elite_exome --extract ${out}/elite_exome.prune.in --chr 1-22 --allow-extra-chr --ibc --out ${out}/elite_exome

# Compute hardy, freq, missingness, genotype counts
plink2 --pfile ${out}/elite_exome --freq --missing --hardy --geno-counts --out ${out}/elite_exome --allow-extra-chr

# filter the data
plink --bfile ${out}/elite_exome --maf 0.001 --hwe 0.001 --geno 0.1 --make-bed --out ${out}/elite_exome.filtered --allow-extra-chr --reference-allele elite.ref.txt
plink2 --bfile ${out}/elite_exome --maf 0.001 --hwe 0.001 --geno 0.1 --make-pgen --out ${out}/elite_exome.filtered --allow-extra-chr --ref-allele elite.ref.txt
