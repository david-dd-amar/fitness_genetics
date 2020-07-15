ml load plink2
ml load plink

SCRIPTS=/oak/stanford/groups/euan/projects/elite/ukbb_exome/scripts/

# set params
# ------------------------------------------------
# Assumptions: these are pgens aligned to the same reference
# -> keeping the ref will avoid any snp flip issues
ukbb=/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome
elite=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_analyze_vcf/elite_exome.filtered

# Ref files
UKREF=/oak/stanford/groups/euan/projects/elite/ukbb_exome/eu_beds/ukbb.pvar.ref
ELREF=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_analyze_vcf/elite.ref.txt

# Missing values
UKGENO=0.01
ELGEINO=0.01

# HWE thresholds
UKHWE=1e-8
ELHWE=1e-6

# MAF thresholds and MAC threhold for the reference
# (included to reduce running time as the ref is very large)
UKMAF=0.01
UKMAC=30
ELMAF=0.01

# Additional MAF data source: compare with gnomAD: remove reference variants
# with a MAF that substantially differs from the ref
gnomADComp=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_filtered_vcfs/comp_freq_with_gnomad
MAX_GNOMAD_DIFF=0.05
MAX_GNOMAD_LOG_OR=0.2

# LD and PCA arguments
LD_PRUNE_R2=0.1
REM_VEXCLUDE=0
PCA_FLAG=""

# Specify population fam for subsetting (done manually)
# POP=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_elite_merged_EU.fam
POP=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_elite_merged_rivaslabEUs.fam
# POP=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_elite_merged_rivaslabEU_intra.fam
# POP=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_elite_merged_rivaslabEU_intra_others.fam
# List of the elite samples
ELSAMPS=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_samples.txt
# Variants to exclude - e.g., intronic or other non exon variants
VEXCLUDE=/oak/stanford/groups/euan/projects/elite/ukbb_exome/variant_anno/non_exon_variants.txt
# Out path
OUT=/oak/stanford/groups/euan/projects/elite/ukbb_exome/eu_beds_v2/
# OUT=/oak/stanford/groups/euan/projects/elite/ukbb_exome/intra_eu_beds_v2/
# OUT=/oak/stanford/groups/euan/projects/elite/ukbb_exome/intra_others_eu_beds_v2/

# additional parameters for GWAS
ELMETA=/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt
UKMETA=/oak/stanford/groups/euan/projects/elite/ukbb_exome/pheno/ukb_pheno_slim.phe

echo "-------------  Parameters for analyzing cohort vs. reference --------------"
echo " Analysis auxiliary scripts are in $SCRIPTS"
echo "  --  cohort path: $elite"
echo "  --  reference path: $ukbb"
echo "  --  remove individuals not in fam file: $POP"
echo "  --  a fam file with our cohort individuals (convinient when splitting data files): $ELSAMPS"
echo "  --  cohort data reference allele per variant: $ELREF"
echo "  --  reference data reference allele per variant: $UKREF"
echo "  --  all results and log files are printed to: $OUT"
echo "------------  parameters for variant filtering  ------------"
echo "  --  cohort missing values filter: --geno $ELGEINO "
echo " 	--  reference missing values filter: --geno $UKGENO "
echo " 	--  cohort hwe filter p-value: $ELHWE"
echo "  --  reference hwe filter p-value: $UKHWE"
echo "  --  for speedup: use minor allele count filter on reference: $UKMAC"
echo "  --  keep variants with MAF>$ELMAF in cohort and MAF>$UKMAF in reference"
echo "  --  compare reference MAFs to nfe MAFs from gnomAD, data path is: $gnomADComp"
echo "  --  remove variants whose reference MAF differs by at least $MAX_GNOMAD_DIFF from gnomAD nfe"
echo " 	--  remove variants whose reference MAF	absolute log OR	with gnomAD nfe is at least $MAX_GNOMAD_LOG_OR"
echo "------------  parameters for LD prune and PCA  ------------"
echo "  --  for PCA and LD prunning, exclude variants from $VEXCLUDE"
echo "  --  LD prunning is done using R2=${LD_PRUNE_R2}"
echo "  --  PCA flag (approx or an empty string): $PCA_FLAG"
echo "  --  Should LD exclude variants from $VEXCLUDE (1-yes, 0-no): $REM_VEXCLUDE"
echo "------------  parameters for GWAS runs  ------------"
echo "  --  cohort metadata file: $ELMETA"
echo "  --  reference metadata file: $UKMETA"

# input vcfs to bed + mac, geno
# -------------------------------------------------------------
# ukbb: extract ref and then make the bed
# slow - run once or use a prev calculated one
# less $ukbb.pvar | cut -f 3,4 > ukbb.pvar.ref
# tail -n +2 ukbb.pvar.ref > ukbb.ref.txt
plink2 --pfile $ukbb vzs --keep $POP --geno $UKGENO --hwe $UKHWE --ref-allele $UKREF --make-bed --mac $UKMAC --out ${OUT}ukbb_exome
# elite
# extract ref file
# slow - run once or use a prev	calculated one
# less $elite | grep -v "#" | cut -f 3,4 > elite.ref.txt
plink2 --pfile $elite --keep $POP --geno $ELGEINO --hwe $ELHWE --make-bed --out ${OUT}elite_exome --allow-extra-chr --rm-dup --ref-allele $ELREF
# 
# merge beds by taking the shared variants
# -------------------------------------------------------------
Rscript ${SCRIPTS}intersect_bims.R ${OUT}ukbb_exome.bim ${OUT}elite_exome.bim $OUT > ${OUT}intersect_bims.log
# to use bmerge we need plink 1.9
ml load plink/1.90b6.16
plink2 --bed ${OUT}ukbb_exome.bed --fam ${OUT}ukbb_exome.fam --bim ${OUT}ukbb_exome.new.bim --extract ${OUT}intersect_variants.txt --ref-allele ${OUT}intersect_variants_ref.txt --out ${OUT}ukbb_exome.intersect --make-bed
plink --bed ${OUT}elite_exome.bed --fam ${OUT}elite_exome.fam --bim ${OUT}elite_exome.new.bim --bmerge ${OUT}ukbb_exome.intersect --allow-extra-chr --out ${OUT}ukbb_elite_merged.tmp
plink2 --bfile ${OUT}ukbb_elite_merged.tmp --extract ${OUT}intersect_variants.txt --ref-allele ${OUT}intersect_variants_ref.txt --out ${OUT}ukbb_elite_merged --make-bed --allow-extra-chr
rm ${OUT}ukbb_elite_merged.tmp*

# subset using freq filters
plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged --freq --keep $ELSAMPS --hardy --geno-counts --out ${OUT}ukbb_elite_merged.elite
plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged --freq --remove $ELSAMPS --hardy --geno-counts --out ${OUT}ukbb_elite_merged.ukbb
Rscript ${SCRIPTS}extract_freq_vars.R ${OUT}ukbb_elite_merged.ukbb.afreq ${OUT}ukbb_elite_merged.elite.afreq $UKMAF $ELMAF > ${OUT}maf_filtered_variants.txt 
plink2 --bfile ${OUT}ukbb_elite_merged --ref-allele ${OUT}intersect_variants_ref.txt --extract ${OUT}maf_filtered_variants.txt --make-bed --out ${OUT}ukbb_elite_merged.filtered

# add analysis of ref vs. gnomad
Rscript ${SCRIPTS}compare_maf_to_gnomad.R $gnomADComp ${OUT}ukbb_elite_merged.ukbb.afreq $MAX_GNOMAD_DIFF $MAX_GNOMAD_LOG_OR > ${OUT}gnomad_maf_excl_variants.txt
plink2 --bfile ${OUT}ukbb_elite_merged.filtered --ref-allele ${OUT}intersect_variants_ref.txt --exclude ${OUT}gnomad_maf_excl_variants.txt --make-bed --out ${OUT}ukbb_elite_merged.filtered
rm *~

# take chr1-22, exclude VEXCLUDE, run LD prunning and PCA
if [ $REM_VEXCLUDE -eq 1 ]
then
	plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged.filtered --exclude $VEXCLUDE --indep-pairwise 500 10 $LD_PRUNE_R2 --out ${OUT}ukbb_elite_merged.filtered
fi
if [ $REM_VEXCLUDE -eq 0 ]
then
        plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged.filtered --indep-pairwise 500 10 $LD_PRUNE_R2 --out ${OUT}ukbb_elite_merged.filtered
fi
# exclude MHC for the PCA
plink2 --extract ${OUT}ukbb_elite_merged.filtered.prune.in --chr 1-5,7-22 --bfile ${OUT}ukbb_elite_merged.filtered --pca 20 ${PCA_FLAG} --out ${OUT}ukbb_elite_merged.filtered

# get missingness (for stats later)
plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged.filtered --missing --keep $ELSAMPS --out ${OUT}ukbb_elite_merged.filtered.elite
plink2 --chr 1-22 --bfile ${OUT}ukbb_elite_merged.filtered --missing --remove $ELSAMPS --out ${OUT}ukbb_elite_merged.filtered.ukbb

# Create the phe file for GWAS
Rscript ${SCRIPTS}create_pheno_file.R ${ELMETA} ${UKMETA} ${OUT}ukbb_elite_merged.filtered.eigenvec ${OUT}merged_pheno.phe

# Run GWAS
plink2  --adjust --allow-extra-chr --bfile ${OUT}ukbb_elite_merged.filtered --chr 1-22 --covar ${OUT}merged_pheno.phe --covar-name sex,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize --logistic hide-covar --out ${OUT}naive_gwas_with_exclude --pheno ${OUT}merged_pheno.phe --pheno-name class --exclude ${VEXCLUDE}
plink2  --adjust --allow-extra-chr --bfile ${OUT}ukbb_elite_merged.filtered --chr 1-22 --covar ${OUT}merged_pheno.phe --covar-name sex,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize --logistic hide-covar --out ${OUT}naive_gwas_no_exclude --pheno ${OUT}merged_pheno.phe --pheno-name class
plink2  --adjust --allow-extra-chr --bfile ${OUT}ukbb_elite_merged.filtered --chr 1-22 --covar ${OUT}merged_pheno.phe --covar-name sex,PC1,PC2,PC3,PC4,PC5,age --covar-variance-standardize --logistic hide-covar --out ${OUT}age_adj_gwas_no_exclude --pheno ${OUT}merged_pheno.phe --pheno-name class
plink2  --adjust --allow-extra-chr --bfile ${OUT}ukbb_elite_merged.filtered --chr 1-22 --covar ${OUT}merged_pheno.phe --covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,age --covar-variance-standardize --logistic hide-covar --out ${OUT}age_adj_gwas_no_exclude_10pcs --pheno ${OUT}merged_pheno.phe --pheno-name class
#plink2  --adjust --allow-extra-chr --bfile ${OUT}ukbb_elite_merged.filtered --chr 1-22 --covar ${OUT}merged_pheno.phe --covar-name sex,PC1,PC2,PC3,PC4,PC5,age --covar-variance-standardize --logistic --out age_adj_gwas_no_exclude --pheno ${OUT}merged_pheno.phe --pheno-name class 

## Run GWAS vs. PCs 1-10
#for pc in {1..10}
#do
#plink2 --chr 1-5,7-22 --bfile ukbb_elite_merged.filtered --pheno-name PC${pc} --pheno ukbb_elite_merged.filtered.eigenvec --linear --adjust --out gwas_pc${pc}
#done

