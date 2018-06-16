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
read_plink_table<-function(path){
  y = read.delim(path,stringsAsFactors = F,header=F)
  y = t(apply(y,1,function(x){
                          x = gsub(x,pattern="^\\s+",replacement = "")
                          strsplit(x,split="\\s+")[[1]]
                        }))
  colnames(y) = y[1,]
  rownames(y) = y[,2]
  return(y[-1,])
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

# Define input parameters
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
snp_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"

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

# create bed file
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



# compare to known sex from the metadata
metadata_sex = sample_metadata_raw$Sex_input_data
names(metadata_sex) = apply(sample_metadata_raw[,1:2],1,paste,collapse="_")
rownames(sample_metadata_raw) = names(metadata_sex)
imputed_sex = read.delim(paste(job_dir,"impute_sex.sexcheck",sep=''),
                         stringsAsFactors = F,header=F)
imputed_sex = t(apply(imputed_sex,1,
                      function(x){
                        x = gsub(x,pattern="^\\s+",replacement = "")
                        strsplit(x,split="\\s+")[[1]]
                      }))
dim(imputed_sex)
colnames(imputed_sex) = imputed_sex[1,]
rownames(imputed_sex) = imputed_sex[,2]
inds = intersect(names(metadata_sex),rownames(imputed_sex))
x1=imputed_sex[inds,4];x2=metadata_sex[inds]
sex_errs = inds[((x1=="1"&x2=="F")|(x1=="2"&x2=="M")) & !is.na(x2)]
sample_metadata_raw[sex_errs,]$Cohort
sex_errs[1:10]
sample_metadata_raw[sex_errs[1:6],]
# look at call rates of the samples with error in sex imputation
quantile(sample_data[sex_errs,]$Call.Rate)
quantile(sample_data$Call.Rate)

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
# Look at the features of the resulting dataset
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

# look at the results
missinigness_report = read_plink_table(paste(job_dir,"maf_filter_missing.imiss",sep=''))
call_rates_after_filters = 1-as.numeric(missinigness_report[,6])
quantile(call_rates_after_filters)
table(call_rates_after_filters<0.97)
low_cr_samples = missinigness_report[call_rates_after_filters<0.98,2]
sample_metadata_raw[low_cr_samples,]$Cohort
intersect(sex_errs,low_cr_samples)

# Exclude samples with low call rate

# Separate chrs 1-22 from other snps

# Infer sex and compare to our covariates

# Recalculate call rates

# Impute missing values

# Run gwas (multiclass, logistic)














