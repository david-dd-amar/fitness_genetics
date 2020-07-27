# This script handles our entire GWAS analysis flow
# It creates a directory with all input, sh, log, err, and output files

####################################################################################################
####################################################################################################
####################################################################################################

# Define analysis parameters for the Illumina reports
# autosomal_chrs = T
snp_min_clustersep_thr = 0.3
snp_min_call_rate = 0.95
snp_min_het_ex = -0.4
snp_max_het_ex = 0.4
# min_maf = 0.001 # For later use
run_loacally = F

# Define analysis parameters for the Illumina reports
initial_subj_min_call_rate = 0.95
final_subj_min_call_rate = 0.98
snp_het_p = 1e-8

script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Define out dir
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_all_batches/"
# set the job's directory
try({system(paste("mkdir",job_dir),wait = T)})
setwd(job_dir)

# 1. MEGG: recalled alone (MEGG: MEGA global)
input_bfile1 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/global_b37_no_reclustering/PLINK_020520_0243/global_b37_no_reclustering"
snp_report_file1 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/global_b37_no_reclustering/SNP_Table.txt"
sample_report_file1 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/global_b37_no_reclustering/Samples_Table.txt"

# 2. MEGC (MEGA Consortium) - reclustred
input_bfile2 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/consortium_b37_with_reclustering/PLINK_100520_1057/consortium_b37_with_reclustering"
snp_report_file2 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/consortium_b37_with_reclustering/SNP_Table.txt"
sample_report_file2 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/consortium_b37_with_reclustering/Samples_Table.txt"

# Additional input
bad_snps_file = "/oak/stanford/groups/euan/projects/fitness_genetics/bad_mega_snps.txt"
build_fa = "/oak/stanford/groups/mrivas/public_data/genomes/hg19/hg19.fa"

library(data.table)
####################################################################################################
####################################################################################################
####################################################################################################
# We need to merge the two MEGA sub datasets that were called separately
# Read the reports and compare
snp_data1 = fread(snp_report_file1,data.table = F,stringsAsFactors = F,check.names = T)
rownames(snp_data1) = snp_data1$Name
sample_data1 = read.delim(sample_report_file1,stringsAsFactors = F)
rownames(sample_data1) = sample_data1$Sample.ID

snp_data2 = fread(snp_report_file2,data.table = F,stringsAsFactors = F,check.names = T)
rownames(snp_data2) = snp_data2$Name
sample_data2 = read.delim(sample_report_file2,stringsAsFactors = F)
rownames(sample_data2) = sample_data2$Sample.ID

# some preprocessing of metadata
analyze_snp_report_get_snps_to_exclude<-function(snp_data,snp_min_call_rate,
    snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
    clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep"){
  snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
  snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
    snp_data[[clust_sep_col_name]] < snp_min_clustersep_thr |
    snp_data$Het.Excess < snp_min_het_ex |
    snp_data$Het.Excess > snp_max_het_ex
  return(snp_data$Name[snps_to_exclude])
}
excl1 = analyze_snp_report_get_snps_to_exclude(snp_data1,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep")
excl2 = analyze_snp_report_get_snps_to_exclude(snp_data2,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep")
excl3 = read.table(bad_snps_file,stringsAsFactors = F)[,1]

print("Initial QC using genomestudios reports")
print(paste("excluded from first file:",length(excl1)))
print(paste("excluded from second file:",length(excl2)))
print(paste("intersect:",length(intersect(excl1,excl2))))

# Final snp set for subsequent analyses
snp_set = intersect(rownames(snp_data1),rownames(snp_data2))
snp_set = setdiff(snp_set,excl1)
snp_set = setdiff(snp_set,excl2)
snp_set = setdiff(snp_set,excl3)
print(paste("using genomestudio data, the number of surviving snps:",length(snp_set)))
write.table(t(t(snp_set)),paste(job_dir,"snp_set.txt",sep=""), row.names = F,col.names = F,quote = F)
snp_file = paste(job_dir,"snp_set.txt",sep="")

# Subject QC: missingness
# File 1
err_path = paste(job_dir,"subj_qc_bfile1.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile1,
                 "--missing --het --make-bed",
                 "--extract", snp_file,
                 "--out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "subj_qc_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# File 2
err_path = paste(job_dir,"subj_qc_bfile2.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile2,
                 "--missing --het --make-bed",
                 "--extract", snp_file,
                 "--out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "subj_qc_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# Analyze the results
# Subject call rates:
miss1 = read.table(paste(job_dir,"bfile1.imiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.imiss",sep=""),stringsAsFactors = F,header=T)
crs1 = 1-miss1$F_MISS;crs2 = 1-miss2$F_MISS
het1 = read.table(paste(job_dir,"bfile1.het",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.het",sep=""),stringsAsFactors = F,header=T)

# Remove problematic samples
to_rem2 = miss2[crs2<initial_subj_min_call_rate | het2$F< -0.25 | het2$F > 0.25 ,1:2]
to_rem1 = miss1[crs1<initial_subj_min_call_rate | het1$F< -0.25 | het1$F > 0.25 ,1:2]
remove_subjects_using_plink(paste(job_dir,"bfile2",sep=""),
                            to_rem2,job_dir,"file2_initial_subj_qc","bfile2",
                            batch_script_func=get_sh_default_prefix)
remove_subjects_using_plink(paste(job_dir,"bfile1",sep=""),
                            to_rem1,job_dir,"file1_initial_subj_qc","bfile1",
                            batch_script_func=get_sh_default_prefix)
wait_for_job(60)
print("After initial qc, datas sizes are:")
print(paste("number of samples, file 1:",length(readLines(paste(job_dir,"bfile1.fam",sep="")))))
print(paste("number of snps, file 1:",length(readLines(paste(job_dir,"bfile1.bim",sep="")))))
ids1 = read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F)[,2]
ids2 = read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of samples, file 2:",length(readLines(paste(job_dir,"bfile2.fam",sep="")))))
print(paste("number of snps, file 2:",length(readLines(paste(job_dir,"bfile2.bim",sep="")))))

# SNP QC: Get missingness, het, call rates etc for each file
# File 1
err_path = paste(job_dir,"snp_c_bfile1.err",sep="")
log_path = paste(job_dir,"snp_c_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--freq 	--hardy --missing",
                 "--out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "snp_c_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# File 2
err_path = paste(job_dir,"snp_c_bfile2.err",sep="")
log_path = paste(job_dir,"snp_c_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
                 "--freq 	--hardy --missing",
                 "--out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "snp_c_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# Analyze the SNP QC results
miss1 = read.table(paste(job_dir,"bfile1.lmiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.lmiss",sep=""),stringsAsFactors = F,header=T)
frq1 = read.table(paste(job_dir,"bfile1.frq",sep=""),stringsAsFactors = F,header=T)
frq2 = read.table(paste(job_dir,"bfile2.frq",sep=""),stringsAsFactors = F,header=T)
het1 = read.table(paste(job_dir,"bfile1.hwe",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.hwe",sep=""),stringsAsFactors = F,header=T)

to_rem_file1 = miss1$SNP[1-miss1$F_MISS<snp_min_call_rate]
to_rem_file1 = union(to_rem_file1,het1$SNP[het1$P < snp_het_p])
to_rem_file2 = miss2$SNP[1-miss2$F_MISS<snp_min_call_rate]
to_rem_file2 = union(to_rem_file2,het2$SNP[het2$P < snp_het_p])

# Some numbers
sum(1-miss1$F_MISS<snp_min_call_rate, na.rm=T)
sum(het1$P < snp_het_p, na.rm=T)
print(paste("Second QC using plink, removed from file 1:",length(to_rem_file1)))
print(paste("Second QC using plink, removed from file 2:",length(to_rem_file2)))
print(paste("Second QC, overlap between the two files:",length(intersect(to_rem_file2,to_rem_file1))))

to_rem_files_union = union(to_rem_file1,to_rem_file2)
bim1 = fread(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F,data.table = F)
snps_to_keep = setdiff(bim1[,2],to_rem_files_union)
print(paste("remaining number of snps after second qc:",length(snps_to_keep)))

# Reduce the datasets before the merge
extract_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile1",
    batch_script_func=get_sh_default_prefix)
extract_snps_using_plink(paste(job_dir,"bfile2",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile2",
                         batch_script_func=get_sh_default_prefix)
wait_for_job()

print("After initial qc, data sizes are:")
print(paste("number of samples, file 1:",length(readLines(paste(job_dir,"bfile1.fam",sep="")))))
print(paste("number of snps, file 1:",length(readLines(paste(job_dir,"bfile1.bim",sep="")))))
ids1 = read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F)[,2]
ids2 = read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of samples, file 2:",length(readLines(paste(job_dir,"bfile2.fam",sep="")))))
print(paste("number of snps, file 2:",length(readLines(paste(job_dir,"bfile2.bim",sep="")))))

# Transform to vcf
# In combination with a FASTA file, --ref-from-fa sets REF alleles when it can be done unambiguously. 
# (Note that this is never possible for deletions and some insertions.)
# File 1
err_path = paste(job_dir,"vcf_bfile1.err",sep="")
log_path = paste(job_dir,"vcf_bfile1.log",sep="")
curr_cmd = c("ml load plink2",paste("plink2 --bfile",paste(job_dir,"bfile1",sep=''),
                 "--export vcf bgz id-paste=iid",
                 "--fa",build_fa,"--ref-from-fa force",
                 "--out",paste(job_dir,"vcfs/bfile1",sep='')))
curr_sh_file = "vcf_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# File 2
err_path = paste(job_dir,"vcf_bfile2.err",sep="")
log_path = paste(job_dir,"vcf_bfile2.log",sep="")
curr_cmd = c("ml load plink2",paste("plink2 --bfile",paste(job_dir,"bfile2",sep=''),
                 "--export vcf bgz id-paste=iid",
                 "--fa",build_fa,"--ref-from-fa force",
                 "--out",paste(job_dir,"vcfs/bfile2",sep='')))
curr_sh_file = "vcf_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

####################################################################################################
####################################################################################################
####################################################################################################
# Merge using PLINK
err_path = paste(job_dir,"merge_plink.err",sep="")
log_path = paste(job_dir,"merge_plink.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--bmerge",paste(job_dir,"bfile2",sep=''),
                 "--make-bed --out",paste(job_dir,"merged_mega_data",sep=''))
curr_sh_file = "merge_plink.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)
print("After merge, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data.bim",sep="")))))

# Check flipscan
# create a phe file
fam1 = as.matrix(read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F))
fam2 = as.matrix(read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F))
fam1 = cbind(fam1,rep("1",nrow(fam1)))
fam2 = cbind(fam2,rep("2",nrow(fam2)))
rownames(fam1)=NULL;colnames(fam1)=NULL
rownames(fam2)=NULL;colnames(fam2)=NULL
phe = rbind(fam1,fam2)
phe = phe[,c(1:2,7)]
colnames(phe) = c("FID","IID","mega")
write.table(phe,file=paste(job_dir,"flipscan_mega_type.txt",sep=""),
            sep=" ",quote=F,row.names = F,col.names = T)
# Run flipscan
err_path = paste(job_dir,"mega_flipscan.err",sep="")
log_path = paste(job_dir,"mega_flipscan.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--flip-scan --allow-no-sex",
                 "--pheno",paste(job_dir,"flipscan_mega_type.txt",sep=""),
                 "--pheno-name mega",
                 "--out",paste(job_dir,"merged_mega_data",sep=''))
curr_sh_file = "mega_flipscan.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(120)

flipscan_res = readLines(paste(job_dir,"merged_mega_data.flipscan",sep=""))
arrs = strsplit(flipscan_res[-1],split="\\s+")
names(arrs) = sapply(arrs,function(x)x[3])
table(sapply(arrs,length))
flipscan_failures = sapply(arrs,length) > 11
flipscan_failures = sapply(arrs[flipscan_failures],function(x)x[3])
# sapply(arrs[flipscan_failures],function(x)x[2])
bim = fread(paste(job_dir,"merged_mega_data.bim",sep=""),stringsAsFactors = F,data.table = F)
snps_to_keep = setdiff(bim[,2],flipscan_failures)
print(paste("Flipscan check, number of variants to remove:",length(flipscan_failures)))
extract_snps_using_plink(paste(job_dir,"merged_mega_data",sep=""),snps_to_keep,job_dir,
                         "final_snps_to_keep_after_flipscan","merged_mega_data",
                         batch_script_func=get_sh_default_prefix)
wait_for_job(60)

print("After flipscan, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data.bim",sep="")))))

####################################################################################################
####################################################################################################
####################################################################################################
# LD, Rel/Kin, PCA, metadata load
err_path = paste(job_dir,"ld_report.err",sep="")
log_path = paste(job_dir,"ld_report.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--indep-pairwise 50 10",0.1,
                 "--chr 1-22",
                 "--out",paste(job_dir,"merged_mega_data",sep=""))
curr_sh_file = "ld_report.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# PCA and King estimates
err_path = paste(job_dir,"pca_king.err",sep="")
log_path = paste(job_dir,"pca_king.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--extract", paste(job_dir,"merged_mega_data.prune.in",sep=''),
                 "--pca",
                 "--make-king-table --king-table-filter 0.177",
                 "--out",paste(job_dir,"merged_mega_data",sep=""))
curr_sh_file = "pca_king.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)


# Analyze the results
# 1. PCA for EU selection: see the "pca_pop_struct_analyses" notebook
# (done manually) - make sure that we have "#IID" as the column name of 
# /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe

# 2. Select rel to exclude 
simple_vc<-function(edges){
  s = c()
  while(length(c(edges))>2){
    topdeg_v = names(sort(table(c(edges)),decreasing=T))[1]
    to_rem = apply(edges==topdeg_v,1,any)
    s = c(s,topdeg_v)
    edges = edges[!to_rem,]
  }
  if(length(edges)>0){
    s = c(s,edges[1])
  }
  return(s)
}
rel_res = fread(
  "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_all_batches/merged_mega_data.kin0",
   stringsAsFactors=F,data.table=F)
rel_edges = as.matrix(rel_res[,c("IID1","IID2")])
rel_to_remove = simple_vc(rel_edges)
rel_to_remove = c("#IID",rel_to_remove)
write.table(t(t(rel_to_remove)),file="/oak/stanford/groups/euan/projects/fitness_genetics/pheno/rel_to_remove.phe",
            row.names=F,col.names=F,sep="\t",quote=F)

# 3. Sex checks
# This is a problematic process, so we use it as a general test and to impute
# missing sex values. We do not exclude samples without additional manual inspection.
# Sex checks were run on the raw data separately using:
# file=global_b37_no_reclustering
# OR
# file=consortium_b37_with_reclustering
# 
# plink --bfile ${file} --chr 23-24 --make-bed --out sexd_unsplit
# plink --bfile sexd_unsplit --split-x b37 no-fail --make-bed --out sexd_split
# plink --bfile sexd_split --indep-pairphase 20000 2000 0.5 --chr 23-24 --out sexd_split
# plink --bfile sexd_split --extract sexd_split.prune.in --make-bed --out sexd_split_pruned
# plink --bfile sexd_split_pruned --check-sex --out ${file}
# Thus for each a raw input bfile we have a .sexcheck file we can read - these are considered
# when creating the master phenotypic table


####################################################################################################
####################################################################################################
####################################################################################################

# Set the fid to zero - this makes --keep and --remove work with a single
# #IID column.
famd = fread(paste(job_dir,"merged_mega_data.fam",sep=''),stringsAsFactors = F,data.table = F)
famd[,1] = 0
write.table(famd,file=paste(job_dir,"merged_mega_data.fam",sep=''),
            sep = " ",row.names = F,col.names = F,quote=F)

# Rerun LD and PCA using EUs
err_path = paste(job_dir,"eu_ld_report.err",sep="")
log_path = paste(job_dir,"eu_ld_report.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--indep-pairwise 50 10",0.1,
                 "--chr 1-22",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "eu_ld_report.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# PCA
err_path = paste(job_dir,"eu_pca.err",sep="")
log_path = paste(job_dir,"eu_pca.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--extract", paste(job_dir,"merged_mega_data.eu.prune.in",sep=''),
                 "--pca",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "eu_pca.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# Repeat the LD and PCA, this time exclude the EU PCA outliers
# from /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe
err_path = paste(job_dir,"eu_ld_report2.err",sep="")
log_path = paste(job_dir,"eu_ld_report2.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--indep-pairwise 50 10",0.1,
                 "--chr 1-22",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu2",sep=""))
curr_sh_file = "eu_ld_report2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# PCA
err_path = paste(job_dir,"eu_pca2.err",sep="")
log_path = paste(job_dir,"eu_pca2.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--extract", paste(job_dir,"merged_mega_data.eu2.prune.in",sep=''),
                 "--pca",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu2",sep=""))
curr_sh_file = "eu_pca2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)
# At this point we have all the data we need:
# Pheno - EUs, PCs, EU PCs, sex checks, failed samples, relatedness results
# Filtered merged data
# Next steps (in different scripts)
#   Create a master phe file
#   Run GWAS
#   Run imputation
#   Analyze the data























####################################################################################################
####################################################################################################
####################################################################################################
# TBD - not sure we need the code below
####################################################################################################
####################################################################################################
####################################################################################################

# Run the GWASs
covs_file = paste(curr_dir,"all_covs_and_pheno.phe",sep="")
for(pcs in 1:5){
  # Vs GP
  curr_pcs = paste(paste("PC",1:pcs,sep=""),collapse = ",")
  cov_line = paste("--covar-name sex,age,",curr_pcs,sep="")
  curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                   "--logistic hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name elite_vs_gp",
                   "--covar",covs_file,
                   "--maf 0.01",
                   cov_line,
                   "--allow-no-sex --adjust",
                   "--threads",8,
                   "--out",paste(curr_dir,"elite_vs_gp_gwas_res_PCs",pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_dir,paste("elite_vs_gp_gwas_res_PCs",pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
  
  curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                   "--logistic hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name cooper_vs_gp",
                   "--covar",covs_file,
                   "--maf 0.01",
                   cov_line,
                   "--allow-no-sex --adjust",
                   "--threads",8,
                   "--out",paste(curr_dir,"cooper_vs_gp_gwas_res_PCs",pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_dir,paste("cooper_vs_gp_gwas_res_PCs",pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
  
  # VO2
  curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                   "--linear hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name VO2max..ml.kg.min.",
                   "--covar",covs_file,
                   "--maf 0.01",
                   cov_line,
                   "--allow-no-sex --adjust",
                   "--threads",8,
                   "--out",paste(curr_dir,"elite_vo2_ml_kg_min_PCs",pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_dir,paste("elite_vo2_ml_kg_min_PCs",pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
  
  curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                   "--linear hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name VO2max..l.",
                   "--covar",covs_file,
                   "--maf 0.01",
                   cov_line,
                   "--allow-no-sex --adjust",
                   "--threads",8,
                   "--out",paste(curr_dir,"elite_vo2max_l_PCs",pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_dir,paste("elite_vo2max_l_PCs",pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
  
  # Add Cooper's run test
  curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                   "--linear hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name Treadmill.time",
                   "--covar",covs_file,
                   "--maf 0.01",
                   cov_line,
                   "--allow-no-sex --adjust",
                   "--threads",8,
                   "--out",paste(curr_dir,"cooper_treadmill_time_PCs",pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_dir,paste("cooper_treadmill_time_PCs",pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
  
}

# Print in fuma-compatible format
out_files = list.files(curr_dir)
out_files = out_files[grepl("logistic$",out_files)]
out_files = out_files[!grepl("fuma",out_files)]
for(ff in out_files){
  print(ff)
  res = fread(paste(curr_dir,ff,sep=""),stringsAsFactors = F,data.table = F)
  print(paste("num variants with p<1e-06:",sum(as.numeric(res[,"P"]<1e-6),na.rm = T)))
  print(paste("num variants with p<5e-08:",sum(as.numeric(res[,"P"]<5e-8),na.rm = T)))
  fuma_res = res[,c("CHR","BP","OR","P")]
  colnames(fuma_res) = c("chromosome","position","OR","P-value")
  write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
              sep=" ")
}
# For linear
out_files = list.files(curr_dir)
out_files = out_files[grepl("linear$",out_files)]
out_files = out_files[!grepl("fuma",out_files)]
for(ff in out_files){
  print(ff)
  res = read.table(paste(curr_dir,ff,sep=""),stringsAsFactors = F,header=T)
  print(paste("num variants with p<1e-06:",sum(as.numeric(res[,"P"]<1e-6),na.rm = T)))
  print(paste("num variants with p<5e-08:",sum(as.numeric(res[,"P"]<5e-8),na.rm = T)))
  fuma_res = res[,c("CHR","BP","BETA","P")]
  colnames(fuma_res) = c("chromosome","position","beta","P-value")
  write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
              sep=" ")
}

# locally: merge elite vo2 and cooper max
setwd("~/Desktop/elite/fuma_results/input_files/")
pcs = 5
library(data.table);library(metap)
elite_res = fread("./fuma_elite_vo2_ml_kg_min_PCs5.assoc.linear",
                  data.table = F,stringsAsFactors = F)
elite_res = elite_res[!is.na(elite_res$`P-value`),]
elite_res = elite_res[!(elite_res[,1]==0 & elite_res[,2]==0),]
rownames(elite_res) = paste(elite_res[,1],elite_res[,2])
cooper_res = fread("./fuma_cooper_treadmill_time_PCs5.assoc.linear",
                   data.table = F,stringsAsFactors = F)
cooper_res = cooper_res[!is.na(cooper_res$`P-value`),]
cooper_res = cooper_res[!(cooper_res[,1]==0 & cooper_res[,2]==0),]
rownames(cooper_res) = paste(cooper_res[,1],cooper_res[,2])
shared_loci = intersect(rownames(cooper_res),rownames(elite_res))
length(shared_loci)
elite_res = elite_res[shared_loci,]
cooper_res = cooper_res[shared_loci,]
ps1 = elite_res$`P-value`
ps2 = cooper_res$`P-value`
merged_ps = apply(cbind(ps1,ps2),1,metap::minimump)
merged_ps = sapply(merged_ps,function(x)x$p)
table(merged_ps<1e-05)
table(p.adjust(merged_ps,method="fdr")<0.5)


####################################################################################################
####################################################################################################
####################################################################################################
# Transform the dataset into HRC-based or 1000G-based data
# Run the check_bim analysis

curr_dir = paste(job_dir,"1000g/",sep="")
system(paste("mkdir",curr_dir))
setwd(curr_dir)
curr_bfile = "merged_mega_data_autosomal"
check_if_bim_is_sorted(paste(paste(job_dir,curr_bfile,".bim",sep='')))
err_path = paste(curr_dir,"run_check_bim.err",sep="")
log_path = paste(curr_dir,"run_check_bim.log",sep="")
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",curr_dir))
# For 1000G-based analysis
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(job_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(job_dir,curr_bfile,".frq",sep=''),
                 "-1000g -p EUR -t 0.3 -r ",
                 "/home/users/davidama/apps/check_bim/1000GP_Phase3_combined.legend")
curr_sh_file = "run_check_bim.sh"
# print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
#               get_sh_prefix_bigmem(err_path,log_path,mem_size = 256000,time="6:00:00"),curr_cmd)
# Try without bigmem
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
wait_for_job()
wait_for_job()
system(paste("less ",curr_dir,"Run-plink.sh | grep TEMP > ",curr_dir,"Run-plink_1000g.sh",sep=""))
run_sh_lines = readLines(paste(curr_dir,"Run-plink_1000g.sh",sep=""))
run_sh_lines[1] = gsub(paste("plink --bfile",curr_bfile),paste("plink --bfile ",job_dir,curr_bfile,sep=""),run_sh_lines[1])
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "")
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

# For HRC-based analysis
curr_dir = paste(job_dir,"hrc/",sep="")
system(paste("mkdir",curr_dir))
setwd(curr_dir)
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",curr_dir))
err_path = paste(curr_dir,"run_check_bim2.err",sep="")
log_path = paste(curr_dir,"run_check_bim2.log",sep="")
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(job_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(job_dir,curr_bfile,".frq",sep=''),
                 "-hrc -p EU -t 0.3 -r",
                 "/home/users/davidama/apps/check_bim/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
curr_sh_file = "run_check_bim2.sh"
# print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
#               get_sh_prefix_bigmem(err_path,log_path,mem_size = 256000,time="6:00:00"),curr_cmd)
# Try without bigmem
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
wait_for_job()
system(paste("less ",curr_dir,"Run-plink.sh | grep TEMP > ",curr_dir,"Run-plink_hrc.sh",sep=""))
run_sh_lines = readLines(paste(curr_dir,"Run-plink_hrc.sh",sep=""))
run_sh_lines[1] = gsub(paste("plink --bfile",curr_bfile),paste("plink --bfile ",job_dir,curr_bfile,sep=""),run_sh_lines[1])
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "")
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

# To download and install the tools on the cluster
# 1. Check bim:
# wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# unzip HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# mkdir check_bim
# mv HRC* check_bim/
# mv LICENSE.txt check_bim/
# cd check_bim
# wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
# gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
#
# 2. Strand analysis:
# mkdir ~/apps/wrayner_strand/
# cd ~/apps/wrayner_strand
# wget http://www.well.ox.ac.uk/~wrayner/strand/update_build.sh
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Use PCA and relatedness to remove samples and then rerun PCA
# covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),stringsAsFactors = F)
# rownames(covariate_matrix) = covariate_matrix$IID
# altsamp_id = sample_metadata_raw$alt_sample_id
# names(altsamp_id) = rownames(sample_metadata_raw)
# is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = rownames(sample_metadata_raw)
# 
# set.seed(123)
# d = covariate_matrix
# pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
# rownames(pc_x) = rownames(d)
# # kmeans_res = kmeans(pc_x,4)$cluster
# # Using hierarchical
# dd = dist(pc_x,method="manhattan")
# h = hclust(dd,method = "single")
# kmeans_res = run_hclust(pc_x,12,dd,h)
# 
# table(kmeans_res)
# table(kmeans_res,d$Cohort)
# table(kmeans_res,is_jap[names(kmeans_res)])
# 
# to_rem = rep(F,nrow(d))
# for(j in 1:20){
#   x = d[,paste("PC",j,sep="")]
#   x = (x-mean(x))/sd(x)
#   print(sum(abs(x)>6))
#   to_rem[abs(x)>6] = T
# }
# table(to_rem)
# 
# # Select subjects from the largest cluster
# clustable = table(kmeans_res)
# selected_cluster = names(which(clustable == max(clustable)))
# selected_subjects = names(which(kmeans_res == selected_cluster))
# 
# print(paste("PCA and clustering analysis, largest cluster size is:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# selected_subjects = setdiff(selected_subjects,rownames(d)[to_rem])
# 
# print(paste("PC outliers removed, number of remaining subjects:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# library("igraph",lib.loc = "~/R/packages")
# rl_data = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf.genome",sep=""),
#                      header=T,stringsAsFactors = F)
# rl_edges = as.matrix(rl_data[,c("IID1","IID2")])
# rl_g = igraph::graph_from_edgelist(rl_edges,directed = F)
# rl_clusters = clusters(rl_g)[[1]]
# rl_subjects_to_remove = c()
# for(cl in unique(rl_clusters)){
#   curr_subjects = names(rl_clusters)[rl_clusters==cl]
#   rl_subjects_to_remove = c(rl_subjects_to_remove,curr_subjects[-1])
# }
# print(paste("Relatedness analysis, number of subjects to remove:",length(rl_subjects_to_remove)))
# print(paste("Intersection with largest cluster",
#             length(intersect(rl_subjects_to_remove,selected_subjects))))
# selected_subjects = setdiff(selected_subjects,rl_subjects_to_remove)
# print(paste("PCA and clustering analysis, largest cluster size is:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# curr_fam = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf.fam",sep=""),stringsAsFactors = F)
# curr_fam = curr_fam[!is.element(curr_fam[,2],set=selected_subjects),]
# remove_subjects_using_plink(paste(job_dir,"merged_mega_data_autosomal_after_maf",sep=""),
#                             curr_fam,
#                             job_dir,"_pca_and_rl_subj_qc","merged_mega_data_autosomal_after_maf_after_pca",
#                             batch_script_func=get_sh_default_prefix)
# wait_for_job()
# print("After PCA and Rl analysis, data sizes are:")
# print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.fam",sep="")))))
# print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.bim",sep="")))))
# 
# # Prune and run PCA and relatedness
# analysis_name = "our_data_ld_prune2"
# err_path = paste(job_dir,analysis_name,"_ld_report2.err",sep="")
# log_path = paste(job_dir,analysis_name,"_ld_report2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--indep-pairwise 250 10",0.1,
#                  "--out",paste(job_dir,analysis_name,"_after_pca",sep=""))
# curr_sh_file = paste(analysis_name,"_ld_report2.sh",sep="")
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# print(paste("Prune before PCA, num of variants is:",
#             length(readLines(paste(job_dir,analysis_name,"_after_pca.prune.in",sep="")))))
# err_path = paste(job_dir,"final_data_pca2.err",sep="")
# log_path = paste(job_dir,"final_data_pca2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--extract", paste(job_dir,analysis_name,"_after_pca.prune.in",sep=""),
#                  "--pca --out",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''))
# curr_sh_file = "final_data_pca2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# # rerun frq
# err_path = paste(job_dir,"final_data_frq2.err",sep="")
# log_path = paste(job_dir,"final_data_frq2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--freq --out",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''))
# curr_sh_file = "final_data_frq2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# 
# # Create the new covariate matrix file
# new_pca_res = read_pca_res(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.eigenvec",sep=""))
# covariate_matrix = covariate_matrix[rownames(new_pca_res),]
# covariate_matrix[,paste("PC",1:20,sep="")] = new_pca_res
# write.table(covariate_matrix,file=
#               paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.txt",sep=''),
#             sep="\t",quote=F,row.names = F)
# write.table(covariate_matrix,file=
#               paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.phe",sep=''),
#             sep=" ",quote=F,row.names = F)
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # output all gwas results into a new dir with input files for fuma interpretation
# 
# res1 = read.table(paste(job_dir,"cooper_treadmill_times.Treadmill.time.glm.linear.adjusted",sep=""),
#                   stringsAsFactors = F)
# res2 = read.table(paste(job_dir,"elite_treadmill_times.vo2.glm.linear.adjusted",sep=""),
#                   stringsAsFactors = F)
# x1 = res1[,3];names(x1)=res1[,2]
# x2 = res2[,3];names(x2)=res2[,2]
# cor(x1,x2[names(x1)])
# fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
# new_ps = apply(cbind(x1,x2[names(x1)]),1,fishersMethod)
# 
# curr_bfile = paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep='')
# create_fuma_files_for_fir(job_dir,
#                           paste(curr_bfile,".bim",sep=""),
#                           paste(curr_bfile,".frq",sep=""),p = 1,maf = 0.01,
#                           snps_to_exclude_from_results=NULL)
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # QC steps: compare our datasets with other datasets
# 
# # for cohort specific mafs see _cohort_freq as regex
# 
# # 1. MAFs: ours vs. the ukbb
# our_mafs = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.frq",sep=""),
#                       header=T,stringsAsFactors = F)
# ukbb_mafs = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq",
#                        header=T,stringsAsFactors = F)
# rownames(our_mafs) = our_mafs$SNP
# rownames(ukbb_mafs) = ukbb_mafs$SNP
# curr_snps = intersect(our_mafs$SNP,ukbb_mafs$SNP)
# x1 = our_mafs[curr_snps,"MAF"]
# x2 = ukbb_mafs[curr_snps,"MAF"]
# cor(x1,x2,method="spearman") # Spearman correlation: 0.997
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 441961    146 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 442034     73 
# 
# # 2. MAFs: ours vs. the ukbb after 1000g as ref
# our_mafs = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.frq",sep=""),
#                       header=T,stringsAsFactors = F)
# ukbb_mafs = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.frq",
#                        header=T,stringsAsFactors = F)
# rownames(our_mafs) = our_mafs$SNP
# rownames(ukbb_mafs) = ukbb_mafs$SNP
# curr_snps = intersect(our_mafs$SNP,ukbb_mafs$SNP)
# x1 = our_mafs[curr_snps,"MAF"]
# x2 = ukbb_mafs[curr_snps,"MAF"]
# cor(x1,x2,method="spearman") # Spearman correlation: 0.997
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 437340     85 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 437395     30
# y1 = our_mafs[curr_snps,c("A1","A2")]
# y1_flipped = t(apply(y1,1,flip_snp_info))
# y2 = ukbb_mafs[curr_snps,c("A1","A2")]
# # table(y1[,1]==y2[,1])
# # FALSE   TRUE 
# # 58274 379151 
# # table(y1[,2]==y2[,1])
# # FALSE   TRUE 
# # 432706   4719 
# # table(y1_flipped[,1]==y2[,1])
# # FALSE   TRUE 
# # 382415  55010
# strand_flipped_inds = y1_flipped[,1]==y2[,1]
# maf_flipped_inds = y1[,2]==y2[,1]
# # quantile(x2[maf_flipped_inds])
# # 0%      25%      50%      75%     100% 
# # 0.001459 0.170950 0.481500 0.494000 0.500000 
# # quantile(x2[strand_flipped_inds])
# # 0%       25%       50%       75%      100% 
# # 0.0009767 0.1382000 0.2552000 0.3751000 0.5000000 
# # quantile(x2)
# # 0%       25%       50%       75%      100% 
# # 0.0007508 0.0513200 0.1662000 0.3174000 0.5000000 
# 
# # Check the JHU stuff
# jhu_snps = our_mafs[grepl(our_mafs$SNP,pattern = "JHU"),]
# our_bim = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.bim",sep=""),
#                       header=F,stringsAsFactors = F)
# rownames(our_bim) = our_bim[,2]
# ukbb_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.bim",
#                        header=F,stringsAsFactors = F)
# inds1 = is.element(ukbb_bim[,4],set = our_bim[jhu_snps$SNP,4])
# ukbb_bim = ukbb_bim[inds1,]
# jhu_bim = our_bim[jhu_snps$SNP,]
# rownames(ukbb_bim) = ukbb_bim[,2]
# ids1 = apply(jhu_bim[,c(1,4)],1,paste,collapse=";")
# ids2 = apply(ukbb_bim[,c(1,4)],1,paste,collapse=";")
# rownames(ukbb_bim) = ids2
# jhu_bim = jhu_bim[ids1!="0;0",]
# ids1 = ids1[ids1!="0;0"]
# rownames(jhu_bim) = ids1
# ids = intersect(ids1,ids2)
# jhu_bim = jhu_bim[ids,]
# ukbb_bim = ukbb_bim[ids,]
# x1 = our_mafs[jhu_bim[,2],"MAF"]
# x2 = ukbb_mafs[ukbb_bim[,2],"MAF"]
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 142586   1977 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 143707    856 


####################################################################################################
####################################################################################################
####################################################################################################
# # Run ld reports (useful for ldscore analysis)
# # File 1
# err_path = paste(job_dir,"ld_bfile1.err",sep="")
# log_path = paste(job_dir,"ld_bfile1.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
#                  "--r2 --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0.1 --freq",
#                  "--out",paste(job_dir,"bfile1",sep=''))
# curr_sh_file = "ld_bfile1.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# # File 2
# err_path = paste(job_dir,"ld_bfile2.err",sep="")
# log_path = paste(job_dir,"ld_bfile2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
#                  "--r2 --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0.1 --freq",
#                  "--out",paste(job_dir,"bfile2",sep=''))
# curr_sh_file = "ld_bfile2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# 
# # Analyze the results
# bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
# bim2 = read.table(paste(job_dir,"bfile2.bim",sep=""),stringsAsFactors = F)
# rownames(bim1) = bim1[,2];rownames(bim2) = bim2[,2]
# 
# ld1 = read.table(paste(job_dir,"bfile1.ld",sep=""),stringsAsFactors = F,header=T)
# ld1_snps = table(c(ld1[,6],ld1[,3]))
# 
# res = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/gwas/gwas_three_groups_linear.ExerciseGroup.glm.linear.adjusted",stringsAsFactors = F)
# res_top_snps = res[res[,3]<1e-50,2]
# ld1_snps[res_top_snps]

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# Locally, should be commented out before running as a batch
setwd("/Users/David/Desktop/elite/july2019_analysis/")
d = read.delim("integrated_sample_metadata_and_covariates.txt",stringsAsFactors = F)
rownames(d) = d$IID
d2 = read.delim("../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
table(d$Cohort)
d$Cohort[d$Cohort=="2"]="ELITE"
d$Cohort[d$Cohort=="1"] = "Cooper"
d$is_jap = is_jap[rownames(d)]

# Cluster by the first two PCs
set.seed(123)
pc_x = as.matrix(d[,paste("PC",1:2,sep="")])
rownames(pc_x) = rownames(d)

to_rem = rep(F,nrow(d))
for(j in 1:10){
  x = d[,paste("PC",j,sep="")]
  x = (x-mean(x))/sd(x)
  print(sum(abs(x)>8))
  to_rem[abs(x)>8] = T
}
table(to_rem)
# pc_x = pc_x[!to_rem,]

# Code for automatic clustering
# pc_x_kmeans = kmeans(pc_x,2)
# kmeans_res = pc_x_kmeans$cluster
# # Using hierarchical
# dd = dist(pc_x,method="euclidean")
# h = hclust(dd,method = "complete")
# kmeans_res = run_hclust(pc_x,5,dd,h)
# kmeans_res[kmeans_res!=1] = 0
# 
# kmeans_res[rownames(d)[to_rem]] = 10
# kmeans_res = kmeans_res[rownames(d)]
# table(kmeans_res)
# table(kmeans_res,d$Cohort)
# table(kmeans_res,is_jap[names(kmeans_res)])
# 
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#   d$PC2[inds],kmeans_res[inds],kmeans_res[inds],
#   main = "PCA+Clustering",xlab="PC1",ylab="PC2",lwd=2,cex.axis=1.4,cex.lab=1.4,
#   xlim = c(-0.02,0.02),ylim = c(-0.05,0.1))
# legend(x="topright",names(res[[1]]),fill = res[[1]],cex=1.3,ncol = 4)

library(ggplot2)

inds = !is.na(d$Cohort)
# inds = 1:nrow(d)
ggplot(d[inds,],aes(x=PC1, y=PC2,shape=Cohort,color=Cohort)) + 
  geom_point(size=2) + ggtitle("PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))
ggplot(d[inds,],aes(x=PC3, y=PC4,shape=Cohort,color=Cohort)) + 
  geom_point(size=2) + ggtitle("PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))
ggplot(d[inds,],aes(x=PC5, y=PC6,shape=Cohort,color=Cohort)) + 
  geom_point(size=2) + ggtitle("PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))
ggplot(d[inds,],aes(x=PC1, y=PC2,shape=Cohort,color=is_jap)) + 
  geom_point(size=2) + ggtitle("PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))

# manual_clustering = d$PC1 < 0 & d$PC2 < 0.02 # Before July 2019
manual_clustering = d$PC1 < 0 & d$PC2 < 0.02 # After July 2019
d$manual_clustering = manual_clustering
ggplot(d[inds,],aes(x=PC1, y=PC2,shape=Cohort,color=manual_clustering)) + 
  geom_point(size=2) + ggtitle("PCs 1 and 2") + 
  theme(plot.title = element_text(hjust = 0.5))

names(manual_clustering) = rownames(d)
table(manual_clustering,to_rem)
table(manual_clustering)
table(manual_clustering,d$Cohort)
table(manual_clustering,is_jap[names(manual_clustering)])
manual_clustering[to_rem] = F
table(manual_clustering)
table(manual_clustering,d$Cohort)

save(manual_clustering,file="manual_clustering.RData")
write.table(t(t(manual_clustering)),file="manual_clustering.txt",sep="\t",
            quote=F)

# # Check number of clusters in PCA plot
# wss <- sapply(1:10,
#               function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
# wss <- sapply(1:20,function(k){tot_wss_hluct(k,h,pc_x)})
# 
# plot(1:length(wss), wss,
#      type="b", pch = 19, frame = FALSE,cex.lab=1.5,lwd=2,cex.axis=1.4,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")
# 
# pc_ps = c()
# for(j in 1:20){
#   # pc_ps[j] = compute_pc_vs_discrete_variable_association_p(d[,paste("PC",j,sep="")],
#   #                                                          d$Cohort)
#   pc_ps[j] = compute_pc_vs_binary_variable_association_p(d[,paste("PC",j,sep="")],
#                                                          d$Cohort)
# }
# p.adjust(pc_ps)

source("~/Desktop/repos/fitness_genetics/R/gwas_flow_helper_functions.R")
pcax = as.matrix(d[,paste("PC",1:40,sep="")])
pcax = pcax[manual_clustering,]
pc_rocs = c()
for(j in 1:40){
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="Cooper"
  p1 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="genepool"
  p2 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "Cooper" | d[,"Cohort"]=="genepool"
  p3 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  pc_rocs = rbind(pc_rocs,c(p1,p2,p3))
}
pc_rocs > 0.75
pc_ps=c()
for(j in 1:40){
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="Cooper"
  p1 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="genepool"
  p2 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "Cooper" | d[,"Cohort"]=="genepool"
  p3 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  pc_ps = rbind(pc_ps,c(p1,p2,p3))
}
pc_ps < 1e-5
