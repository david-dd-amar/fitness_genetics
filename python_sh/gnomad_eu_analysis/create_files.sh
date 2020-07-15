#!/bin/bash

ml load bcftools

gnomad=/oak/stanford/groups/euan/projects/elite/external_data/gnomad/

# for elite
our=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_filtered_vcfs/
out=/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_filtered_vcfs/gnomad_intersect/

# for ukbb
our=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_analyze_pgen/gnomad_intersect/ukbb_vcfs/
out=/oak/stanford/groups/euan/projects/elite/ukbb_exome/ukbb_analyze_pgen/gnomad_intersect/

RSCRIPT=/oak/stanford/groups/euan/projects/elite/ukbb_exome/scripts/gnomad_eu_analysis/extract_eu_info.R
#out=/oak/stanford/groups/euan/projects/elite/external_data/gnomad/eu_mafs/

for i in {1..22}
do
	mkdir ${out}chr${i}
	bcftools isec -p ${out}chr${i} -n=2 -O z -w1 ${gnomad}gnomad.genomes.r3.0.sites.chr${i}.vcf.bgz ${our}chr${i}.vcf.gz &
done

for i in {1..22}
do
	mv ${out}chr${i}/0000.vcf ${out}chr${i}.vcf
done

for i in {1..22}
do
	Rscript $RSCRIPT ${out}chr${i}.vcf ${out}chr${i}.eu_info.tsv &
done
