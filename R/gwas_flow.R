# This script handles our entire GWAS analysis flow
# It creates a directory with all input, sh, log, err, and output files

####################################################################################################
####################################################################################################
####################################################################################################
######## Some useful functions
get_sh_default_prefix<-function(err="",log=""){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=24:00:00",
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=2",
      "#SBATCH --mem=64000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}
print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}

get_my_jobs<-function(){
  tmp = paste("tmp",abs(rnorm(1)),sep='')
  system(paste("squeue | grep davidama >",tmp),wait = T)
  jobs = readLines(tmp)
  system(paste("rm",tmp))
  return(jobs)
}

get_job_id<-function(x){
  x = strsplit(x,split="\\s+",perl=T)[[1]][-1]
  return(x[1])
}

wait_for_job<-function(jobs_before,waittime=30){
  Sys.sleep(waittime)
  jobs_before = sapply(jobs_before,get_job_id)
  curr_jobs = get_my_jobs()
  curr_jobs = sapply(curr_jobs,get_job_id)
  new_job = setdiff(curr_jobs,jobs_before)
  if(length(new_job)==0){return(NULL)}
  print(paste("new added job is: ",new_job))
  while(is.element(new_job,set=curr_jobs)){
    Sys.sleep(waittime)
    curr_jobs = get_my_jobs()
    curr_jobs = sapply(curr_jobs,get_job_id)
  }
}

correct_dups_in_sample_metadata<-function(x){
  nns = apply(x[,1:2],1,paste,collapse="_")
  dups = names(which(table(nns)>1))
  to_keep = rep(F,length(nns))
  for(i in 1:length(nns)){
    if(!is.element(nns[i],set=dups)){to_keep[i]=T;next}
    curr_inds = which(nns==nns[i])
    num_nas = apply(is.na(x[curr_inds,]),1,sum)
    curr_inds = curr_inds[num_nas==min(num_nas)]
    to_keep[curr_inds[1]]=T
  }
  x = x[to_keep,]
  rownames(x) = nns[to_keep]
  return(x)
}
read_plink_table<-function(path,has_header=T,...){
  y = read.delim(path,stringsAsFactors = F,header=F,...)
  y = t(apply(y,1,function(x){
                          x = gsub(x,pattern="^\\s+",replacement = "")
                          strsplit(x,split="\\s+")[[1]]
                        }))
  rownames(y) = y[,2]
  if(has_header){
    colnames(y) = y[1,]
    return(y[-1,])
  }
  return(y)
}

two_d_plot_visualize_covariate<-function(x1,x2,cov1,cov2=NULL,cuts=5,...){
  if(is.null(cov2)){cov2=cov1}
  if(is.numeric(cov1)){cov1=cut(cov1,breaks = cuts)}
  if(is.numeric(cov2)){cov1=cut(cov2,breaks = cuts)}
  cov1 = as.factor(cov1)
  cov2 = as.factor(cov2)
  cols = rainbow(length(unique(cov1)))
  names(cols) = unique(cov1)
  pchs = 1:length(unique(cov2))
  names(pchs) = unique(cov2)
  plot(x1,x2,col=cols[cov1],pch=pchs[cov2],...)
  return(list(cols,pchs))
}

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

# Define Input parameters
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
snp_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"

# TODO:
# 1. (Later, low pref for now) Adapt the code to handle NULL snp and sample reports.
#    In these cases we do SNP/Sample filtering using PLINK's algorithms
# 3. Before the GWAS: exclude samples with either low quality scores after SNP filters or those
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
wait_for_job(jobs_before,5)
list.files(job_dir)

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
wait_for_job(jobs_before,5)
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
wait_for_job(jobs_before,5)
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
wait_for_job(jobs_before,5)
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
                 "--make-bed --out",paste(job_dir,"maf_filter",sep=''))
curr_sh_file = "maf_filter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
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
wait_for_job(jobs_before,5)
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
wait_for_job(jobs_before,5)
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

####################################################################################################
####################################################################################################
####################################################################################################
# Locally, should be commented out before running as a batch
setwd("/Users/David/Desktop/elite/analysis/")
d = read.delim("integrated_sample_metadata_and_covariates.txt")

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

inds = d$Cohort !="genepool"
res = two_d_plot_visualize_covariate(d$PC2[inds],d$PC3[inds],d$Cohort[inds],d$Cohort[inds])
legend(x="topleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d$PC2[inds],d$PC3[inds],
                                     d$Cohort[inds],d$Cohort[inds],
                                     xlim = c(0,0.01),
                                     ylim = c(-0.01,0.01))


####################################################################################################
####################################################################################################
####################################################################################################
# From here: GWAS
# Exclude samples with low call rate and error in sex inference

# Impute missing values (?)

# Create phe file for GWAS

# Run gwas (multiclass, logistic)














