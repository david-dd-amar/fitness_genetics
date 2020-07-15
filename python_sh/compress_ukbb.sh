#!/bin/bash

#SBATCH --mem 10G
#SBATCH --time=24:00:00
#SBATCH -c 8
##SBATCH --account=mrivas
#SBATCH -p euan,mrivas,normal,owners
#SBATCH --output=compress_ukbb.out
#SBATCH --error=compress_ukbb.err
##SBATCH -p nih_s10

ml load gatk
ml load bcftools

bgzip --threads 8 ukbb_exome.vcf
