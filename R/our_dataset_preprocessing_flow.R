# This script handles our entire GWAS analysis flow
# It creates a directory with all input, sh, log, err, and output files

####################################################################################################
####################################################################################################
####################################################################################################

# Define analysis parameters
autosomal_chrs = T
snp_min_clustersep_thr = 0.5
snp_min_call_rate = 0.95
snp_min_het_ex = -0.3
snp_max_het_ex = 0.2
min_maf = 0.005
run_loacally = F
num_pca_clusters = 3
analysis_cohorts = "elite" # If null then use all cohorts

# Define Input parameters
# Original PLINK files
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
# ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
# map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
# New fwd strand files (August 2018)
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/"
# Alternative for elite only: August 2018
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/our_prepro/"


ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.ped"
map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.map"

snp_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)
setwd(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# read reports and metadata
snp_data = read.delim(snp_report_file)
sample_data = read.delim(sample_report_file,stringsAsFactors = F)
rownames(sample_data) = sample_data$Sample.ID
sample_metadata_raw = read.delim(sample_metadata,stringsAsFactors = F)
sample_metadata_raw = correct_dups_in_sample_metadata(sample_metadata_raw)

# some preprocessing of metadata
snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
     snp_data$Multi.EthnicGlobal_D1.bpm.Cluster.Sep < snp_min_clustersep_thr |
     snp_data$Het.Excess < snp_min_het_ex |
     snp_data$Het.Excess > snp_max_het_ex
if(autosomal_chrs){
  snps_to_exclude = snps_to_exclude & snp_data_autosomal_rows
}
if(!autosomal_chrs){
  snps_to_exclude = snps_to_exclude & !snp_data_autosomal_rows
}
table(snps_to_exclude)
write.table(t(t(as.character(snp_data$Name[snps_to_exclude]))),
            file=paste(job_dir,"snps_to_exclude.txt",sep=''),
            row.names = F,col.names = F,quote = F)

####################################################################################################
####################################################################################################
####################################################################################################
# set the job's directory
system(paste("mkdir",job_dir),wait = T)
system(paste("cp",ped_file,paste(job_dir,"raw.ped",sep='')),wait = T)
system(paste("cp",map_file,paste(job_dir,"raw.map",sep='')),wait = T)
list.files(job_dir)

# Create bed file
jobs_before = get_my_jobs()
err_path = paste(job_dir,"raw_to_bed.err",sep="")
log_path = paste(job_dir,"raw_to_bed.log",sep="")
curr_cmd = paste("plink --file",paste(job_dir,"raw",sep=''),
                 "--make-bed --out",paste(job_dir,"raw",sep=''))
curr_sh_file = "raw_bed.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

####################################################################################################
####################################################################################################
####################################################################################################
# For analysis of specific cohorts
if(!is.null(analysis_cohorts)){
  sample_metadata_raw$Cohort = tolower(sample_metadata_raw$Cohort)
  analysis_cohorts = tolower(analysis_cohorts)
  curr_samples = rownames(sample_metadata_raw)[is.element(sample_metadata_raw$Cohort,
                                                          set=analysis_cohorts)]
  
  fam_samples = read.table(paste(job_dir,"raw.fam",sep=""),stringsAsFactors = F)
  iid2fid = fam_samples[,1]
  names(iid2fid) = fam_samples[,2]
  curr_samples = intersect(curr_samples,fam_samples[,2])
  curr_samples = cbind(iid2fid[curr_samples],curr_samples)
  write.table(curr_samples,file=paste(job_dir,"curr_samples.txt",sep=""),
              row.names = F,quote = F,col.names = F)
  
  err_path = paste(job_dir,"keep_cohort_subjects.err",sep="")
  log_path = paste(job_dir,"keep_cohort_subjects.log",sep="")
  curr_cmd = paste("plink --bfile",paste(job_dir,"raw",sep=''),
                   "--keep",paste(job_dir,"curr_samples.txt",sep=''),
                   "--make-bed --out",paste(job_dir,"raw",sep=''))
  curr_sh_file = "keep_cohort_subjects.sh"
  print_sh_file(paste(job_dir,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
  wait_for_job()
  system(paste("rm ",job_dir,"raw.bed~",sep=""))
  system(paste("rm ",job_dir,"raw.bim~",sep=""))
  system(paste("rm ",job_dir,"raw.fam~",sep=""))
}

####################################################################################################
####################################################################################################
####################################################################################################
# impute sex
err_path = paste(job_dir,"impute_sex.err",sep="")
log_path = paste(job_dir,"impute_sex.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"raw",sep=''),
                 "--check-sex .3 .9 --chr 23",
                 "--out",paste(job_dir,"impute_sex",sep=''))
curr_sh_file = "impute_sex.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# get simple imputation by looking at call rates in Y chromosome
Y_snps = grepl("^Y$",snp_data$Chr,ignore.case = T)
write.table(t(t(as.character(snp_data$Name[Y_snps]))),
            file=paste(job_dir,"Y_snps.txt",sep=''),
            row.names = F,col.names = F,quote = F)
jobs_before = get_my_jobs()
err_path = paste(job_dir,"Y_snps_analysis.err",sep="")
log_path = paste(job_dir,"Y_snps_analysis.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"raw",sep=''),
                 "--extract", paste(job_dir,"Y_snps.txt",sep=""),
                 "--missing --out",paste(job_dir,"Y_snps_analysis",sep=''))
curr_sh_file = "Y_snps_analysis.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################
# keep the desired chrs (default: 0-22)
# also, exclude low qc score snps 
# in theory we can use the ped itself for inference and filtering
# however, here we assume we have the snp report that we can use
# above the list of excluded snps was printed to a file in the job dir
err_path = paste(job_dir,"chr_filter.err",sep="")
log_path = paste(job_dir,"chr_filter.log",sep="")
chr_filter = "--chr 0-22"
if(!autosomal_chrs){
  chr_filter = "--no-chr 0-22"
}
curr_cmd = paste("plink --bfile",paste(job_dir,"raw",sep=''),
                 chr_filter,
                 "--exclude",paste(job_dir,"snps_to_exclude.txt",sep=''),
                 "--make-bed --out",paste(job_dir,"chr_filter",sep=''))
curr_sh_file = "chr_filter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# quick comparison of the data using TOP strand vs + strand
# plus_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/raw.bim")
# top_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/raw.bim")
# dim(plus_bim)
# dim(top_bim)
# comparisons = list()
# for(j in 1:ncol(plus_bim)){
#   comparisons[[j]] = table(plus_bim[,j]==top_bim[,j])  
# }
####################################################################################################
####################################################################################################
####################################################################################################
# Exclude low maf snps.
# Also add freq, pca, and misingness analyses
if(is.null(analysis_cohorts)){
  maf_snps_to_exclude = snp_data$Name[snp_data$Minor.Freq<min_maf]
  write.table(t(t(as.character(snp_data$Name[maf_snps_to_exclude]))),
              file=paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
              row.names = F,col.names = F,quote = F)
  err_path = paste(job_dir,"maf_filter.err",sep="")
  log_path = paste(job_dir,"maf_filter.log",sep="")
  curr_cmd = paste("plink --bfile",paste(job_dir,"chr_filter",sep=''),
                   "--exclude",paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
                   "--maf",min_maf,
                   "--pca --freq --missing",
                   "--make-bed --out",paste(job_dir,"maf_filter_data",sep=''))
  curr_sh_file = "maf_filter.sh"
  print_sh_file(paste(job_dir,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
  wait_for_job()
}

if(!is.null(analysis_cohorts)){
  err_path = paste(job_dir,"maf_filter.err",sep="")
  log_path = paste(job_dir,"maf_filter.log",sep="")
  curr_cmd = paste("plink --bfile",paste(job_dir,"chr_filter",sep=''),
                   "--maf 0.01", # for elite, consider changing for other cohorts
                   "--pca --freq --missing",
                   "--make-bed --out",paste(job_dir,"maf_filter_data",sep=''))
  curr_sh_file = "maf_filter.sh"
  print_sh_file(paste(job_dir,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
  wait_for_job()
}

####################################################################################################
####################################################################################################
####################################################################################################
# Here we analyze the results of the jobs above
# Needed files:
#   Y-based sex analysis: Y_snps_analysis.imiss
#   Plink's sex analysis: impute_sex.sexcheck
# Analyze the sex analysis results
Y_missing_report = read_plink_table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/Y_snps_analysis.imiss")
Y_inferred_sex = Y_missing_report[,6] == "nan"
# compare to known sex from the metadata
metadata_sex = sample_metadata_raw$Sex_input_data
names(metadata_sex) = apply(sample_metadata_raw[,1:2],1,paste,collapse="_")
rownames(sample_metadata_raw) = names(metadata_sex)
imputed_sex = read_plink_table(paste(job_dir,"impute_sex.sexcheck",sep=''))
inds = intersect(names(metadata_sex),rownames(imputed_sex))
x1=imputed_sex[inds,4];x2=metadata_sex[inds];x3 = Y_inferred_sex[inds]
sex_errs1 = inds[((x1=="1"&x2=="F")|(x1=="2"&x2=="M")) & !is.na(x2)]
sex_errs2 = inds[((!x3 & x2=="F")|(x3 & x2=="M")) & !is.na(x2)]
sex_errs = union(sex_errs1,sex_errs2)
# Check the ids and compare to the exome data
# Create a report with the sex errors
m = sample_metadata_raw[sex_errs,]$Cohort
# look at call rates of the samples with error in sex imputation
m = cbind(m,sample_data[sex_errs,]$Call.Rate)
m = cbind(x2[sex_errs],x1[sex_errs],x3[sex_errs],m)
m = cbind(sample_metadata_raw[sex_errs,"Sample_ID"],m)
rownames(m) = sex_errs
colnames(m) = c("Sample_ID","Sex_metadata","Sex_plink_impute(1_is_male)","Sex_y_inf(is_female)","Cohort","Raw_call_rate")
write.table(m,file=paste(job_dir,"sex_impute_analysis_report.txt",sep=''),sep="\t",quote=F)

# Missigness report after quality and maf filtering
# look at the results, compare to Illumina's 
missinigness_report = read_plink_table(paste(job_dir,"maf_filter_data.imiss",sep=''))
call_rates_after_filters = 1-as.numeric(missinigness_report[,6])
names(call_rates_after_filters) = rownames(missinigness_report)
low_cr_samples = missinigness_report[call_rates_after_filters<0.98,2]
raw_call_rates = sample_data[rownames(missinigness_report),"Call.Rate"]

# Report low call rate samples
m = sample_metadata_raw[low_cr_samples,]$Cohort
# look at call rates of the samples with error in sex imputation
m = cbind(m,sample_data[low_cr_samples,]$Call.Rate)
m = cbind(m,call_rates_after_filters[low_cr_samples])
m = cbind(sample_metadata_raw[low_cr_samples,"Sample_ID"],m)
colnames(m) = c("Sample_ID","Cohort","Raw_call_rate","Call_rate_qc_maf_filters")
write.table(m,file=paste(job_dir,"missigness_analysis_report.txt",sep=''),sep="\t",quote=F)

# Analyze the PCA results
pca_res = read_plink_table(paste(job_dir,"maf_filter_data.eigenvec",sep=''),F)
pca_res = pca_res[,-c(1:2)]
colnames(pca_res) = paste("PC",1:ncol(pca_res),sep='')

# some sanity checks before printing the output report
all(rownames(pca_res) == names(call_rates_after_filters))
all(rownames(pca_res) == names(imputed_sex))
all(rownames(pca_res) == names(raw_call_rates))

# define the samples to exclude for subsequent analysis, get the bed, bgen, and covariate files
fam_samples = read.table(paste(job_dir,"raw.fam",sep=""),stringsAsFactors = F)
subjects_for_analysis = intersect(rownames(sample_metadata_raw),fam_samples[,2])
subjects_for_analysis = setdiff(subjects_for_analysis,low_cr_samples)
subjects_for_analysis = setdiff(subjects_for_analysis,sex_errs)

# put all covariates in one table
covariate_matrix = cbind(sample_metadata_raw[rownames(pca_res),c(6:8,11,12:17,20:23)],
                         raw_call_rates,call_rates_after_filters,imputed_sex,
                         pca_res)
# note that the plink data may have additional samples not in the samples
# we wish to analyze in this specific project.
# We therefore have to make sure we proceed only with subjects in the metadata
# file.
covariate_matrix = covariate_matrix[subjects_for_analysis,]
length(intersect(sex_errs,rownames(covariate_matrix)))
write.table(covariate_matrix,file=
              paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),
            sep="\t",quote=F)

####################################################################################################
####################################################################################################
####################################################################################################

write.table(covariate_matrix[,c("FID","IID")],file=
              paste(job_dir,"subjects_for_analysis.txt",sep=''),
            sep="\t",quote=F,row.names = F)
err_path = paste(job_dir,"exclude_failed_subjects.err",sep="")
log_path = paste(job_dir,"exclude_failed_subjects.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"maf_filter_data",sep=''),
                 "--keep",paste(job_dir,"subjects_for_analysis.txt",sep=''),
                 "--pca --freq --missing --make-bed --out",paste(job_dir,"final_dataset_for_analysis",sep=''))
curr_sh_file = "exclude_failed_subjects.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################
# Create a freq file for each cohort
covariate_matrix = read.table(paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),sep="\t")
table(covariate_matrix$Cohort)
for(cc in unique(covariate_matrix$Cohort)){
  inds = covariate_matrix$Cohort == cc
  m = covariate_matrix[inds,c("FID","IID")]
  write.table(m,sep=" ",file=paste(job_dir,cc,"_subjects.txt",sep=""),row.names = F,quote = F)
  err_path = paste(job_dir,cc,"_cohort_freq.err",sep="")
  log_path = paste(job_dir,cc,"_cohort_freq.log",sep="")
  curr_cmd = paste("plink --bfile",paste(job_dir,"final_dataset_for_analysis",sep=''),
                   "--keep",paste(job_dir,cc,"_subjects.txt",sep=""),
                   "--freq --out",paste(job_dir,cc,"_cohort_freq",sep=""))
  curr_sh_file = paste(cc,"_cohort_freq.sh",sep="")
  print_sh_file(paste(job_dir,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
}

####################################################################################################
####################################################################################################
####################################################################################################
# Transform the dataset into HRC-based data
# Run the check_bim analysis
setwd(job_dir)
err_path = paste(job_dir,"run_check_bim.err",sep="")
log_path = paste(job_dir,"run_check_bim.log",sep="")
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",job_dir))
curr_cmd = paste("perl", paste(job_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(job_dir,"final_dataset_for_analysis.bim",sep=''),
                 "-f", paste(job_dir,"final_dataset_for_analysis.frq",sep=''),
                 "-hrc -p EUR -r",
                 "/home/users/davidama/apps/check_bim/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
curr_sh_file = "run_check_bim.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# system(paste("mv /home/users/davidama/apps/check_bim/*final_dataset_for_analysis*",job_dir))
# system(paste("mv /home/users/davidama/apps/check_bim/Run-plink.sh",job_dir))
wait_for_job()
system(paste("less ",job_dir,"Run-plink.sh | grep TEMP > ",job_dir,"Run-plink2.sh",sep=""))

err_path = paste(job_dir,"run_check_bim_update.err",sep="")
log_path = paste(job_dir,"run_check_bim_update.log",sep="")
plink_commands = readLines(paste(job_dir,"Run-plink2.sh",sep=""))
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),plink_commands)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

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

# # ####################################################################################################
# # ####################################################################################################
# # ####################################################################################################
# # Locally, should be commented out before running as a batch
# setwd("/Users/David/Desktop/elite/analysis/")
# d = read.delim("integrated_sample_metadata_and_covariates.txt")
# d2 = read.delim("../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
# d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
# samp_id = d2$Sample_ID
# altsamp_id = d2$alt_sample_id
# names(samp_id) = d2_ids
# names(altsamp_id) = d2_ids
# is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
# 
# # Cluster by the first two PCs
# pc_x = as.matrix(d[,c("PC1","PC2")])
# pc_x_5means = kmeans(pc_x,5)
# table(pc_x_5means$cluster)
# table(pc_x_5means$cluster,d$Cohort)
# kmeans_res = pc_x_5means$cluster
# 
# inds = d$Cohort !="genepool"
# inds = 1:nrow(d)
# inds = kmeans_res==2 | kmeans_res==5
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#     d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
#     main = "Cooper and Elite",xlab="PC1",ylab="PC2")
# legend(x="topleft",names(res[[1]]),fill = res[[1]])
# 
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#     d$PC2[inds],kmeans_res[inds],kmeans_res[inds],
#     main = "PCA+Kmeans",xlab="PC1",ylab="PC2")
# legend(x="topleft",names(res[[1]]),fill = res[[1]])
# 
# res = two_d_plot_visualize_covariate(d$PC3[inds],
#     d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
#     main = "Cooper and Elite",xlab="PC3",ylab="PC2")
# legend(x="topleft",names(res[[1]]),fill = res[[1]])
# 
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#     d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds],
#     main = "Cooper and Elite",xlab="PC1",ylab="PC2")
# legend(x="topleft",names(res[[1]]),fill = res[[1]])
# 
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#     d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds],
#     main = "All samples by shipment date",xlab="PC1",ylab="PC2")
# legend(x="topleft",names(res[[1]]),fill = res[[1]])
# 
# curr_is_jap = is_jap[rownames(d)]
# table(curr_is_jap)
# inds = curr_is_jap[inds]
# res = two_d_plot_visualize_covariate(d$PC6[inds],
#     d$PC7[inds],curr_is_jap[inds],curr_is_jap[inds],
#     main = "Is JA sample?",xlab="PC6",ylab="PC7",
#     cex = 1+1*curr_is_jap[inds])
# legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
# 
# two_d_plot_visualize_covariate_ggplot(d$PC1[inds],
#   d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds])
# 
# # Check number of clusters in PCA plot
# wss <- sapply(1:5, 
#               function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
# wss
# plot(1:5, wss,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")

