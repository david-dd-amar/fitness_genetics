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

# Define Input parameters
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
snp_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# TODO:
# 1. (Later, low pref for now) Adapt the code to handle NULL snp and sample reports.
#    In these cases we do SNP/Sample filtering using PLINK's algorithms
# 3. (?) Before the GWAS: exclude samples with either low quality scores after SNP filters or those
#    that failed the sex check.
# 5. Define different GWAS flows and get covariates

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
# impute sex
# select snps to exclude
X_snps = grepl("^X$",snp_data$Chr,ignore.case = T)
X_snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
  snp_data$Multi.EthnicGlobal_D1.bpm.Cluster.Sep < snp_min_clustersep_thr |
  snp_data$Het.Excess < snp_min_het_ex |
  snp_data$Het.Excess > snp_max_het_ex
X_snps_to_exclude = X_snps_to_exclude & X_snps
write.table(t(t(as.character(snp_data$Name[X_snps_to_exclude]))),
            file=paste(job_dir,"X_snps_to_exclude.txt",sep=''),
            row.names = F,col.names = F,quote = F)
# run imputation in plink
jobs_before = get_my_jobs()
err_path = paste(job_dir,"impute_sex.err",sep="")
log_path = paste(job_dir,"impute_sex.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"raw",sep=''),
                 "--exclude X_snps_to_exclude.txt",
                 "--check-sex --out",paste(job_dir,"impute_sex",sep=''))
curr_sh_file = "impute_sex.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()
list.files(job_dir)
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
wait_for_job()
list.files(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# keep the desired chrs (default: 0-22)
# also, exclude low qc score snps 
# in theory we can use the ped itself for inference and filtering
# however, here we assume we have the snp report that we can use
# above the list of excluded snps was printed to a file in the job dir
jobs_before = get_my_jobs()
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
wait_for_job()
list.files(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# Exclude low maf snps
maf_snps_to_exclude = snp_data$Name[snp_data$Minor.Freq<min_maf]
write.table(t(t(as.character(snp_data$Name[maf_snps_to_exclude]))),
            file=paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
            row.names = F,col.names = F,quote = F)
jobs_before = get_my_jobs()
err_path = paste(job_dir,"maf_filter.err",sep="")
log_path = paste(job_dir,"maf_filter.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"chr_filter",sep=''),
                 "--exclude",paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
                 "--maf",min_maf,
                 "--make-bed --out",paste(job_dir,"maf_filter",sep=''))
curr_sh_file = "maf_filter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()
list.files(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# Get missigness results
jobs_before = get_my_jobs()
err_path = paste(job_dir,"maf_filter_missing.err",sep="")
log_path = paste(job_dir,"maf_filter_missing.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"maf_filter",sep=''),
                 "--missing --out",paste(job_dir,"maf_filter_missing",sep=''))
curr_sh_file = "maf_filter_missing.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()
list.files(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# Run PCA
jobs_before = get_my_jobs()
err_path = paste(job_dir,"maf_filter_pca.err",sep="")
log_path = paste(job_dir,"maf_filter_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"maf_filter",sep=''),
                 "--pca --out",paste(job_dir,"maf_filter",sep=''))
curr_sh_file = "maf_filter_pca.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()
list.files(job_dir)

# Run freq
jobs_before = get_my_jobs()
err_path = paste(job_dir,"maf_filter_freq.err",sep="")
log_path = paste(job_dir,"maf_filter_freq.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"maf_filter",sep=''),
                 "--freq --out",paste(job_dir,"maf_filter",sep=''))
curr_sh_file = "maf_filter_freq.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()
list.files(job_dir)


####################################################################################################
####################################################################################################
####################################################################################################
# Here we analyze the results of the jobs above
# Needed files:
#   Y-based sex analysis: Y_snps_analysis.imiss
#   Plink's sex analysis: impute_sex.sexcheck
# Analyze the sex analysis results
Y_missing_report = read_plink_table(paste(job_dir,"Y_snps_analysis.imiss",sep=''))
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
missinigness_report = read_plink_table(paste(job_dir,"maf_filter_missing.imiss",sep=''))
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
pca_res = read_plink_table(paste(job_dir,"maf_filter.eigenvec",sep=''),F)
pca_res = pca_res[,-c(1:2)]
colnames(pca_res) = paste("PC",1:ncol(pca_res),sep='')

# some sanity checks before printing the output report
all(rownames(pca_res) == names(call_rates_after_filters))
all(rownames(pca_res) == names(imputed_sex))
all(rownames(pca_res) == names(raw_call_rates))

# put all covariates in one table
covariate_matrix = cbind(sample_metadata_raw[rownames(pca_res),c(6:8,11,12:17,20:23)],
                         raw_call_rates,call_rates_after_filters,imputed_sex,
                         pca_res)
# note that the plink data may have additional samples not in the samples
# we wish to analyze in this specific project.
# We therefore have to make sure we proceed only with subjects in the metadata
# file.
covariate_matrix = covariate_matrix[intersect(rownames(sample_metadata_raw),
                                              rownames(covariate_matrix)),]
length(intersect(sex_errs,rownames(covariate_matrix)))
write.table(covariate_matrix,file=
              paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),
            sep="\t",quote=F)

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Locally, should be commented out before running as a batch
setwd("/Users/David/Desktop/elite/analysis/")
d = read.delim("integrated_sample_metadata_and_covariates.txt")
d2 = read.delim("../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids

# PCA plots
two_d_plot_visualize_covariate<-function(x1,x2,cov1,cov2=NULL,cuts=5,...){
  if(is.null(cov2)){cov2=cov1}
  if(is.numeric(cov1)){cov1=cut(cov1,breaks = cuts)}
  if(is.numeric(cov2)){cov1=cut(cov2,breaks = cuts)}
  cov1 = as.factor(cov1)
  cov2 = as.factor(cov2)
  cols = rainbow(length(unique(cov1)))
  names(cols) = unique(cov1)
  cols = cols[!is.na(names(cols))]
  pchs = 1:length(unique(cov2))
  names(pchs) = unique(cov2)
  pchs = pchs[!is.na(names(pchs))]
  plot(x1,x2,col=cols[cov1],pch=pchs[cov2],...)
  return(list(cols,pchs))
}

# library(ggplot2)
# two_d_plot_visualize_covariate_ggplot<-function(x1,x2,cov1,cov2=NULL,cuts=5,...){
#   if(is.null(cov2)){cov2=cov1}
#   if(is.numeric(cov1)){cov1=cut(cov1,breaks = cuts)}
#   if(is.numeric(cov2)){cov1=cut(cov2,breaks = cuts)}
#   cov1 = as.factor(cov1)
#   cov2 = as.factor(cov2)
#   cols = rainbow(length(unique(cov1)))
#   names(cols) = unique(cov1)
#   cols = cols[!is.na(names(cols))]
#   pchs = 1:length(unique(cov2))
#   names(pchs) = unique(cov2)
#   pchs = pchs[!is.na(names(pchs))]
#   df = data.frame(x1,x2,col=cols[cov1],pch=pchs[cov2])
#   plot1 = ggplot(aes(x=x1, y=x2, group = cols[cov1]),data = df) + 
#     geom_point(shape=pchs[cov2],alpha=0.2, aes(colour = cols[cov1])) +
#     theme_bw(25)
#   plot(plot1)
#   return(list(cols,pchs))
# }

inds = d$Cohort !="genepool"
inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC1",ylab="PC2")
legend(x="topleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC3",ylab="PC2")
legend(x="topleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds],
    main = "Cooper and Elite",xlab="PC1",ylab="PC2")
legend(x="topleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds],
    main = "All samples by shipment date",xlab="PC1",ylab="PC2")
legend(x="topleft",names(res[[1]]),fill = res[[1]])

curr_is_jap = is_jap[rownames(d)]
table(curr_is_jap)
inds = curr_is_jap[inds]
res = two_d_plot_visualize_covariate(d$PC6[inds],
    d$PC7[inds],curr_is_jap[inds],curr_is_jap[inds],
    main = "Is JA sample?",xlab="PC6",ylab="PC7",
    cex = 1+1*curr_is_jap[inds])
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])


two_d_plot_visualize_covariate_ggplot(d$PC1[inds],
  d$PC2[inds],d$Shipment.date[inds],d$Shipment.date[inds])

####################################################################################################
####################################################################################################
####################################################################################################
# job_dir = "/Users/David/Desktop/elite/analysis/" # for local tests
# run_loacally = T
# From here: prepare data and run GWAS
covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),stringsAsFactors = F)
covariate_matrix[covariate_matrix==""] = NA

# Fill in some info:
# No shipment date into one batch
covariate_matrix$Shipment.date[is.na(covariate_matrix$Shipment.date)] = "uknown"

# read our fam file
fam_info = read_plink_table(paste(job_dir,"maf_filter.fam",sep=""),has_header = F)
iid_to_fid = fam_info[,1]

# Exclude samples with low call rate and error in sex inference
sex_analysis_report = read.delim(file=paste(job_dir,"sex_impute_analysis_report.txt",sep=''))
excluded_samples = union(rownames(sex_analysis_report),rownames(covariate_matrix)[covariate_matrix$call_rates_after_filters<0.98])
remaining_samples = rownames(covariate_matrix)[!is.element(rownames(covariate_matrix),set=excluded_samples)]

# Print excluded into file
m1 = cbind(iid_to_fid[excluded_samples],excluded_samples)
colnames(m1) = c("FID","IID")
write.table(m1,file = paste(job_dir,"gwas_excluded_samples.txt",sep=''),row.names = F,col.names = T,quote=F,sep=" ")

Y_missing_report = read_plink_table(paste(job_dir,"Y_snps_analysis.imiss",sep=''))
Y_inferred_sex = Y_missing_report[,6] == "nan"
Sex = Y_inferred_sex[remaining_samples]
Sex[Sex] = "2"
Sex[Sex=="FALSE"] = "1"
table(Sex)

# Create phe file for GWAS 
pheno_cols = c(
  "Cohort",
  "Shipment.date",
  "Age..at.test.",
  paste("PC",1:10,sep="")
)
pheno_data = cbind(remaining_samples,Sex,covariate_matrix[remaining_samples,pheno_cols])
# Correct the cohort
pheno_data$Cohort[pheno_data$Cohort=="ELITE"] = "1"
pheno_data$Cohort[pheno_data$Cohort=="Cooper"] = "0"
pheno_data$Cohort[pheno_data$Cohort=="genepool"] = "-1"
pheno_data$Cohort = as.numeric(pheno_data$Cohort)

for(j in 4:4){
  pheno_data[[j]] = cov_phe_col_to_plink_numeric_format(pheno_data[[j]])
}
pheno_data[pheno_data$ExerciseGroup==1,]$Age

# write phe file
pheno_data = cbind(as.character(iid_to_fid[remaining_samples]),pheno_data)
colnames(pheno_data) = c("FID","IID","sex","ExerciseGroup","Batch","Age",paste("PC",1:10,sep=""))
write.table(file=paste(job_dir,"three_group_analysis_genepool_controls.phe",sep=''),
            pheno_data,sep=" ",row.names = F,col.names = T,quote=F)

write.table(file=paste(job_dir,"three_group_analysis_sex_update.phe",sep=''),
            pheno_data[,1:3],sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# Clean the data for the analysis


####################################################################################################
####################################################################################################
####################################################################################################
# A set of GWAS runs starts here

# 0. Linear for all groups + age, sex, 5 PCs
pheno_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
covar_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
#jobs_before = get_my_jobs()
err_path = paste(job_dir,"gwas_three_groups_linear.err",sep="")
log_path = paste(job_dir,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# 1. Linear for all groups + sex, 5 PCs
pheno_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
covar_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
#jobs_before = get_my_jobs()
err_path = paste(job_dir,"gwas_three_groups_linear.err",sep="")
log_path = paste(job_dir,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"gwas_three_groups_linear_no_age",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# 3. Logistic Cooper vs. Genepool
table(pheno_data$ExerciseGroup)
sample_inds = pheno_data$ExerciseGroup != "1"
curr_pheno = pheno_data[sample_inds,]
curr_pheno[,4] = as.character(as.numeric(curr_pheno[,4]+1))
pheno_file = paste(job_dir,"cooper_vs_genepool.phe",sep='')
write.table(file=pheno_file,curr_pheno[,c(1:2,4)],sep=" ",row.names = F,col.names = T,quote=F)
covar_file = paste(job_dir,"cooper_vs_genepool_covar.phe",sep='')
write.table(file=covar_file,curr_pheno[,-4],sep=" ",row.names = F,col.names = T,quote=F)
#jobs_before = get_my_jobs()
err_path = paste(job_dir,"cooper_vs_genepool.err",sep="")
log_path = paste(job_dir,"cooper_vs_genepool.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--logistic hide-covar firth-fallback",
                 "--1",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"cooper_vs_genepool",sep=''))
curr_sh_file = "cooper_vs_genepool.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# 1. Logistic with age: elite vs. cooper, 5 PCs
# #jobs_before = get_my_jobs()
# if (!run_loacally){
#   err_path = paste(job_dir,"genepool_controls_simple_linear_wo_age.err",sep="")
#   log_path = paste(job_dir,"genepool_controls_simple_linear_wo_age.log",sep="")
#   curr_cmd = paste(paste(job_dir,"plink2",sep=""),
#                    "--bfile",paste(job_dir,"maf_filter",sep=''),
#                    "--logistic hide-covar firth-fallback",
#                    paste("--pheno",pheno_file),
#                    paste("--pheno-name ExerciseGroup"),
#                    "--allow-no-sex",
#                    "--1",
#                    paste("--covar",covar_file),
#                    "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5,PC6",
#                    "--adjust",
#                    "--out",paste(job_dir,"genepool_controls_simple_linear_wo_age",sep=''))
#   curr_sh_file = "genepool_controls_simple_linear_wo_age.sh"
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
#   #wait_for_job()
#   list.files(job_dir)
#   readLines(err_path)
# }
# if(run_loacally){
#   curr_cmd = paste(paste(job_dir,"plink2",sep=""),
#                    "--bfile",paste(job_dir,"maf_filter",sep=''),
#                    "--logistic hide-covar firth-fallback",
#                    paste("--pheno",pheno_file),
#                    paste("--pheno-name ExerciseGroup"),
#                    "--allow-no-sex",
#                    "--1",
#                    paste("--covar",covar_file),
#                    "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5,PC6",
#                    "--adjust",
#                    "--out",paste(job_dir,"genepool_controls_simple_linear_wo_age",sep=''))
#   system(curr_cmd)
# }
#   
# res_files = list.files(job_dir)
# res_files = res_files[grepl("genepool_controls_simple_linear_wo_age",res_files)]
# res_file = res_files[grepl("adjusted$",res_files) & grepl("logistic",res_files)]
# res = read.delim(paste(job_dir,res_file,sep=''),stringsAsFactors = F)
# res[1:5,]
# hist(res$UNADJ)
# qqplot(-log10(sample(res$UNADJ)[1:10000]),-log10(runif(10000)));abline(0,1)
# res_fuma = from_our_sol_to_fuma_res(res_file,
#                                     paste(job_dir,"maf_filter.bim",sep=""),
#                                     paste(job_dir,"maf_filter.frq",sep=""),maf = 0.01)
# qqplot(-log10(sample(as.numeric(res_fuma[,3]))[1:10000]),-log10(runif(10000)));abline(0,1)
# dim(res_fuma)
# gc()
# write.table(res_fuma,
#             file= paste(job_dir,"genepool_controls_simple_linear_wo_age_fuma_res.txt",sep=""),
#             row.names = F,col.names = T,quote = F,sep=" ")
# 
# 









