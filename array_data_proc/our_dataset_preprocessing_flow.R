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
bad_snps_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/bad_mega_snps.txt"
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

# Define the set of bad snps
bad_snps = readLines(bad_snps_file)
bim_info = fread(paste(job_dir,"merged_mega_data.bim",sep=''),stringsAsFactors = F,data.table = F)
jhu_snps = bim_info[grepl("JHU",bim_info[,2]),2]
poly_snps = bim_info[
  (bim_info[,5] == "A" & bim_info[,6] == "T") | 
    (bim_info[,5] == "T" & bim_info[,6] == "A") |
    (bim_info[,5] == "G" & bim_info[,6] == "C") |
    (bim_info[,5] == "C" & bim_info[,6] == "G") , 2
  ]

all_bad_snps = union(union(bad_snps,jhu_snps),poly_snps)
write.table(t(t(all_bad_snps)),
            file="/oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt",
            row.names = F,col.names = F,quote = F)
write.table(t(t(poly_snps)),
            file="/oak/stanford/groups/euan/projects/fitness_genetics/analysis/poly_snps.txt",
            row.names = F,col.names = F,quote = F)
write.table(t(t(jhu_snps)),
            file="/oak/stanford/groups/euan/projects/fitness_genetics/analysis/jhu_snps.txt",
            row.names = F,col.names = F,quote = F)



####################################################################################################
####################################################################################################
####################################################################################################

# LD, Rel/Kin, PCA, metadata load
err_path = paste(job_dir,"ld_report.err",sep="")
log_path = paste(job_dir,"ld_report.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--exclude /oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt",
                 "--geno 0.01",
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
                 "--exclude /oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt",
                 "--geno 0.01",
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

# Define the PCA outliers along the EU PCs
# use 6 SDs as a threshold
eu_pca_outliers  = c()
eu_pca = fread(paste(job_dir,"merged_mega_data.eu.eigenvec",sep=""),stringsAsFactors = F,data.table = F)
rownames(eu_pca) = eu_pca$IID
eu_pca = eu_pca[,-c(1:2)]
print(head(eu_pca))
for (pc in 1:10){
  pcv = eu_pca[[paste0("PC",pc)]]
  names(pcv) = rownames(eu_pca)
  outls = names(pcv)[(abs(pcv-mean(pcv,na.rm=T)))/sd(pcv,na.rm=T) > 6]
  print(length(outls))
  outls = outls[!is.na(outls)]
  eu_pca_outliers = union(eu_pca_outliers,outls)
}
eu_pca_outliers = c("IID",eu_pca_outliers)
print(length(eu_pca_outliers))
write.table(t(t(eu_pca_outliers)),file="/oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
            row.names=F,col.names=F,sep="\t",quote=F)

# Repeat the LD and PCA, this time exclude the EU PCA outliers
# from /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe
err_path = paste(job_dir,"eu_ld_report2.err",sep="")
log_path = paste(job_dir,"eu_ld_report2.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--exclude /oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt",
                 "--geno 0.01",
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

# Check if we have additional outliers
eu_pca_outliers = readLines("/oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe")[-1]
eu_pca = fread(paste(job_dir,"merged_mega_data.eu2.eigenvec",sep=""),stringsAsFactors = F,data.table = F)
rownames(eu_pca) = eu_pca$IID
eu_pca = eu_pca[,-c(1:2)]
print(head(eu_pca))
sum_new_outliers = 0
for (pc in 1:10){
  pcv = eu_pca[[paste0("PC",pc)]]
  names(pcv) = rownames(eu_pca)
  outls = names(pcv)[(abs(pcv-mean(pcv,na.rm=T)))/sd(pcv,na.rm=T) > 6]
  outls = outls[!is.na(outls)]
  print(length(outls))
  sum_new_outliers = sum_new_outliers + length(outls)
  eu_pca_outliers = union(eu_pca_outliers,outls)
}
eu_pca_outliers = c("IID",eu_pca_outliers)
print(length(eu_pca_outliers))
if(sum_new_outliers >0){
  write.table(t(t(eu_pca_outliers)),file="/oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
              row.names=F,col.names=F,sep="\t",quote=F)
}

# Repeat the LD and PCA, this time exclude the EU PCA outliers
# from /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe
err_path = paste(job_dir,"eu_ld_report3.err",sep="")
log_path = paste(job_dir,"eu_ld_report3.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--maf 0.01",
                 "--exclude /oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt",
                 "--geno 0.01",
                 "--indep-pairwise 50 10",0.1,
                 "--chr 1-22",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu3",sep=""))
curr_sh_file = "eu_ld_report3.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# PCA
err_path = paste(job_dir,"eu_pca3.err",sep="")
log_path = paste(job_dir,"eu_pca3.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--extract", paste(job_dir,"merged_mega_data.eu3.prune.in",sep=''),
                 "--pca",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(job_dir,"merged_mega_data.eu3",sep=""))
curr_sh_file = "eu_pca3.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))


####################################################################################################
####################################################################################################
####################################################################################################
# Transform the dataset into HRC-based or 1000G-based data
# TODO: Create an eu bed file before running

# Run the check_bim analysis

check_bim_script = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/check_bim/HRC-1000G-check-bim-NoReadKey.pl"
check_bim_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/check_bim/"

####################################
# 1000G-based analysis
####################################
curr_dir = paste(job_dir,"checkbim_1000g_eu/",sep="")
system(paste("mkdir -p",curr_dir))
setwd(curr_dir)

# Get the EU bed file
err_path = paste(curr_dir,"eu_bed.err",sep="")
log_path = paste(curr_dir,"eu_bed.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--make-bed",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(curr_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "eu_bed.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
readLines(log_path)

err_path = paste(curr_dir,"addfreq.err",sep="")
log_path = paste(curr_dir,"addfreq.log",sep="")
curr_cmd = paste("plink --bfile",paste(curr_dir,"merged_mega_data.eu",sep=""),
                 "--freq",
                 "--out",paste(curr_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "addfreq.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

curr_bfile = "merged_mega_data.eu"
check_if_bim_is_sorted(paste(paste(curr_dir,curr_bfile,".bim",sep='')))
err_path = paste(curr_dir,"run_check_bim.err",sep="")
log_path = paste(curr_dir,"run_check_bim.log",sep="")
system(paste("cp",check_bim_script,curr_dir))
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(curr_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(curr_dir,curr_bfile,".frq",sep=''),
                 "-1000g -p EUR -t 0.3 -r ",
                 paste0(check_bim_path,"1000GP_Phase3_combined.legend"))
curr_sh_file = "run_check_bim.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

run_sh_lines = readLines(paste(curr_dir,"Run-plink.sh",sep=""))
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

####################################
# For HRC-based analysis
####################################
curr_dir = paste(job_dir,"checkbim_hrc_eu/",sep="")
system(paste("mkdir -p",curr_dir))
setwd(curr_dir)

# Get the EU bed file
err_path = paste(curr_dir,"eu_bed.err",sep="")
log_path = paste(curr_dir,"eu_bed.log",sep="")
curr_cmd = paste("plink2 --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--make-bed",
                 "--keep /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_ids.phe",
                 "--remove /oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe",
                 "--out",paste(curr_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "eu_bed.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
readLines(log_path)

err_path = paste(curr_dir,"addfreq.err",sep="")
log_path = paste(curr_dir,"addfreq.log",sep="")
curr_cmd = paste("plink --bfile",paste(curr_dir,"merged_mega_data.eu",sep=""),
                 "--freq",
                 "--out",paste(curr_dir,"merged_mega_data.eu",sep=""))
curr_sh_file = "addfreq.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

curr_bfile = "merged_mega_data.eu"
#check_if_bim_is_sorted(paste(paste(curr_dir,curr_bfile,".bim",sep='')))
err_path = paste(curr_dir,"run_check_bim.err",sep="")
log_path = paste(curr_dir,"run_check_bim.log",sep="")
system(paste("cp",check_bim_script,curr_dir))
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(curr_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(curr_dir,curr_bfile,".frq",sep=''),
                 "-hrc -p EUR -t 0.3 -r ",
                 paste0(check_bim_path,"HRC.r1-1.GRCh37.wgs.mac5.sites.tab"))
curr_sh_file = "run_check_bim.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
                                                         mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

run_sh_lines = readLines(paste(curr_dir,"Run-plink.sh",sep=""))
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

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
# # Print in fuma-compatible format
# out_files = list.files(curr_dir)
# out_files = out_files[grepl("logistic$",out_files)]
# out_files = out_files[!grepl("fuma",out_files)]
# for(ff in out_files){
#   print(ff)
#   res = fread(paste(curr_dir,ff,sep=""),stringsAsFactors = F,data.table = F)
#   print(paste("num variants with p<1e-06:",sum(as.numeric(res[,"P"]<1e-6),na.rm = T)))
#   print(paste("num variants with p<5e-08:",sum(as.numeric(res[,"P"]<5e-8),na.rm = T)))
#   fuma_res = res[,c("CHR","BP","OR","P")]
#   colnames(fuma_res) = c("chromosome","position","OR","P-value")
#   write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
#               sep=" ")
# }
# # For linear
# out_files = list.files(curr_dir)
# out_files = out_files[grepl("linear$",out_files)]
# out_files = out_files[!grepl("fuma",out_files)]
# for(ff in out_files){
#   print(ff)
#   res = read.table(paste(curr_dir,ff,sep=""),stringsAsFactors = F,header=T)
#   print(paste("num variants with p<1e-06:",sum(as.numeric(res[,"P"]<1e-6),na.rm = T)))
#   print(paste("num variants with p<5e-08:",sum(as.numeric(res[,"P"]<5e-8),na.rm = T)))
#   fuma_res = res[,c("CHR","BP","BETA","P")]
#   colnames(fuma_res) = c("chromosome","position","beta","P-value")
#   write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
#               sep=" ")
# }

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
####################################################################################################
####################################################################################################
####################################################################################################
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

