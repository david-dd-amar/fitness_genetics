# This script is based on https://github.com/mhguo1/TRAPD
# Use python 2.7
# Required packages:
#   optparse
#   operator
#   pybedtools
#
# Run in environment:
# ml load python/2.7.13
# ml load bedtools
# ml load py-pybedtools

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Current assumption: no need to do step 0 of TRAPD:
# Separating multi-allelic and left-aligning indels, annotation, and read depth filter
# see the repo page for details.
# Also, this has to be gzipped
# Taken from here and gzipped:
# vcf_input_path = "/oak/stanford/groups/euan/projects/elite/annotation/elite_n267_joint_snpeff.vcf"
vcf_input_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/elite_n267_joint_snpeff.vcf.gz"
trapd_repo = "/home/users/davidama/apps/trapd/TRAPD/code/"
out_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/"
setwd(out_dir)

# Step 1: create the SNP file
# exmple:
# python2.7 /home/users/davidama/apps/trapd/TRAPD/code/make_snp_file.py --vcffile /oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/elite_n267_joint_snpeff.vcf.gz --outfile /oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/trapd_step1_out.txt --genecolname ANN
cmd = paste(
  "python2.7",
  paste(trapd_repo,"make_snp_file.py",sep=""),
  "--vcffile",vcf_input_path,
  "--outfile",paste(out_dir,"trapd_step1_out.txt",sep=""),
  "--genecolname ANN"
)
system(cmd)
print(length(readLines(paste(out_dir,"trapd_step1_out.txt",sep=""))))

# Step 2: count carriers in the chohort
cmd = paste(
  "python2.7",
  paste(trapd_repo,"count_cases.py",sep=""),
  "--vcffile",vcf_input_path,
  "--minAN 3",
  "--snpfile",paste(out_dir,"trapd_step1_out.txt",sep=""),
  "--outfile",paste(out_dir,"trapd_step2_out.txt",sep="")
)
system(cmd)

# list.files(out_dir)
# system2("ls",args = paste(out_dir,"-l"))

# gnomad link:
# https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
# exac link:
# ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
gnomad_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz"
exac_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/exome_seq/ExAC.r1.sites.vep.vcf.gz"
# Count carriers in ref data
cmd = paste(
  "python2.7",
  paste(trapd_repo,"count_controls.py",sep=""),
  "--vcffile",gnomad_file,
  "--snpfile",paste(out_dir,"trapd_step1_out.txt",sep=""),
  "--database gnomad",
  "--pop FIN,NFE",
  "--outfile",paste(out_dir,"gnomad_counts.txt",sep="")
)
system(cmd)
cmd = paste(
  "python2.7",
  paste(trapd_repo,"count_controls.py",sep=""),
  "--vcffile",exac_file,
  "--snpfile",paste(out_dir,"trapd_step1_out.txt",sep=""),
  "--database exac",
  "--pop FIN,NFE",
  "--outfile",paste(out_dir,"exac_counts.txt",sep="")
)
system(cmd)

# Run TRAPD
# This requires: install.packages("argparse")
# that is, argparse needs to be loaded within the script
cmd = paste(
  "Rscript",
  paste(trapd_repo,"burden.R",sep=""),
  "--casefile",paste(out_dir,"trapd_step2_out.txt",sep=""),
  "--controlfile",paste(out_dir,"exac_counts.txt",sep=""),
  "--casesize 267",
  "--controlsize 36000",
  "--outfile",paste(out_dir,"exac_trapd_output.txt",sep="")
)
system(cmd)

cmd = paste(
  "Rscript",
  paste(trapd_repo,"burden.R",sep=""),
  "--casefile",paste(out_dir,"trapd_step2_out.txt",sep=""),
  "--controlfile",paste(out_dir,"gnomad_counts.txt",sep=""),
  "--casesize 267",
  "--controlsize 22000",
  "--outfile",paste(out_dir,"gnomad_trapd_output.txt",sep="")
)
system(cmd)

# QQ plots
cmd = paste(
  "Rscript",
  paste(trapd_repo,"QQ.R",sep=""),
  "--pvalfile",paste(out_dir,"exac_trapd_output.txt",sep=""),
  "--plotfile",paste(out_dir,"exac_trapd_output.png",sep="")
)
system(cmd)

cmd = paste(
  "Rscript",
  paste(trapd_repo,"QQ.R",sep=""),
  "--pvalfile",paste(out_dir,"gnomad_trapd_output.txt",sep=""),
  "--plotfile",paste(out_dir,"gnomad_trapd_output.png",sep="")
)
system(cmd)

# Read the output in
library("data.table",lib.loc = "~/R/packages")
d = fread(paste(out_dir,"exac_trapd_output.txt",sep=""),
               stringsAsFactors = F,data.table = F)
d = d[d$CASE_TOTAL_AC > 50,]
d2 = fread(paste(out_dir,"gnomad_trapd_output.txt",sep=""),
          stringsAsFactors = F,data.table = F)
d2 = d2[d2$CASE_TOTAL_AC > 50,]

all(d$GENE == d2$GENE)
q1 = p.adjust(d$P_REC)
q2 = p.adjust(d2$P_REC)
cor(q1,q2)
inds = q1<0.001 & q2<0.001
genes = d[inds,1]
genes[grepl(genes,pattern = "myo",ignore.case = T)]
genes[grepl(genes,pattern = "myh",ignore.case = T)]


