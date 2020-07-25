#!/bin/bash

PROGNAME=compare_vcfs
VERSION="1.0"

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	merge vcfs using bcftools then compare using plink.
	  Author: David Amar
	
	Usage: $PROGNAME <vcf1.gz 1><vcf.gz 2><out>
	Notes:
	  - VCFs are merged, taking the variants that have the same REF/ALT
	  - VCF paths should be full, not relative
          - Create a phe file for plink
          - Compare using plink glm
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments"
    usage
    exit 1
fi

vcf1=${1}
vcf2=${2}
OUT=${3}

mkdir -p ${OUT}
cd ${OUT}
ml load bcftools
ml load plink2

# create a phe file
bcftools query -l ${vcf1} | awk '{print $1"\t1"}' > tmp.phe
bcftools query -l ${vcf2} | awk '{print $1"\t2"}' >> tmp.phe
(echo -e "IID\tgroup" ; cat tmp.phe) > samples.phe

bcftools isec -p . -n=2 -O z ${vcf1} ${vcf2}

bcftools index --csi 0000.vcf.gz
bcftools index --csi 0001.vcf.gz
# -m none: only records with identical REF and ALT alleles are compatible
bcftools merge -m none -O z 0000.vcf.gz 0001.vcf.gz -o merged.vcf.gz

# create a pgen file
plink2 --vcf merged.vcf.gz --make-pgen --not-chr 0 --out merged
# run king estimation
plink2 --pfile merged --indep-pairwise 500 10 0.1 --out merged --maf 0.01
plink2 --pfile merged --extract merged.prune.in --make-king-table --king-table-filter 0.4 --out merged
# take the duplicated samples
#cut -f 1 merged.kin0 > tmp.samples
#cut -f 2 merged.kin0 >> tmp.samples
#less tmp.samples | sort | uniq > dup_samples.txt
Rscript comp_vcfs_extract_dups_from_king.R  samples.phe merged.kin0 dup_samples.txt 
rm tmp.samples
rm tmp.phe
# run the gwas
plink2 --pfile merged --glm --pheno-name group --pheno samples.phe --require-pheno --keep dup_samples.txt --out merged --adjust
