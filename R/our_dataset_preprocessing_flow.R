# This script handles our entire GWAS analysis flow
# It creates a directory with all input, sh, log, err, and output files

####################################################################################################
####################################################################################################
####################################################################################################

# Define analysis parameters for the Illumina reports
autosomal_chrs = T
snp_min_clustersep_thr = 0.3
snp_min_call_rate = 0.95
snp_min_het_ex = -0.4
snp_max_het_ex = 0.4
min_maf = 0.001
run_loacally = F
num_pca_clusters = 3
# analysis_cohorts = "elite" # If null then use all cohorts
# analysis_cohorts = NULL # If null then use all cohorts

# Define analysis parameters for the Illumina reports
initial_subj_min_call_rate = 0.95
final_subj_min_call_rate = 0.98
snp_het_p = 1e-4

script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Define Input parameters
# Original PLINK files
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
# ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
# map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
# New fwd strand files (August 2018)
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/"
# Alternative for elite only: August 2018
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/our_prepro/"
# September 2018
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/"
# set the job's directory
try({system(paste("mkdir",job_dir),wait = T)})
setwd(job_dir)

# Each recalling has a set of parameters
ped_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.ped"
map_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.map"
input_bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering"
snp_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
sample_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"

ped_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data.ped"
map_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data.map"
input_bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data"
snp_report_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/SNP_Table.txt"
sample_report_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/Samples_Table.txt"

bad_snps_file = "/scratch/groups/euan/projects/stanford3k/plink/config/remove_snps.txt"

# Our metadata
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
sample_metadata_raw = read.delim(sample_metadata,stringsAsFactors = F)
sample_metadata_raw = correct_dups_in_sample_metadata(sample_metadata_raw)
# As of September 2018 we do not have genepool's metadata: we ignore this until we get it
sample_metadata_raw = sample_metadata_raw[sample_metadata_raw$Cohort!="genepool",]
sample_metadata_raw = sample_metadata_raw[!is.na(sample_metadata_raw[,1]),]
rownames(sample_metadata_raw) = apply(sample_metadata_raw[,1:2],1,paste,collapse="_")

####################################################################################################
####################################################################################################
####################################################################################################
# We need to merge the two MEGA sub datasets that were called separately
# Read the reports and compare
snp_data1 = read.delim(snp_report_file1)
sample_data1 = read.delim(sample_report_file1,stringsAsFactors = F)
rownames(sample_data1) = sample_data1$Sample.ID

snp_data2 = read.delim(snp_report_file2)
sample_data2 = read.delim(sample_report_file2,stringsAsFactors = F)
rownames(sample_data2) = sample_data2$Sample.ID

# Compare some stats
rownames(snp_data1) = snp_data1$Name
rownames(snp_data2) = snp_data2$Name
inds = intersect(rownames(snp_data1),rownames(snp_data2))
snp_samp = sample(inds)[1:200000]
cor(snp_data1[snp_samp,]$Call.Freq,snp_data2[snp_samp,]$Call.Freq,method="spearman")
cor(snp_data1[snp_samp,]$Het.Excess,snp_data2[snp_samp,]$Het.Excess,method="spearman")
snp_samp = snp_samp[!is.na(snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep)]
cor(snp_data1[snp_samp,]$Multi.EthnicGlobal_D1.bpm.Cluster.Sep,snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep,method="spearman")
table(snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep>0.3,snp_data1[snp_samp,]$Multi.EthnicGlobal_D1.bpm.Cluster.Sep>0.3)

# some preprocessing of metadata
analyze_snp_report_get_snps_to_exclude<-function(snp_data,autosomal_chrs,snp_min_call_rate,
                                                 snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
                                                 clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep"){
  snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
  snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
    snp_data[[clust_sep_col_name]] < snp_min_clustersep_thr |
    snp_data$Het.Excess < snp_min_het_ex |
    snp_data$Het.Excess > snp_max_het_ex
  if(autosomal_chrs){
    snps_to_exclude = snps_to_exclude & snp_data_autosomal_rows
  }
  if(!autosomal_chrs){
    snps_to_exclude = snps_to_exclude & !snp_data_autosomal_rows
  }
  return(snp_data$Name[snps_to_exclude])
}
excl1 = analyze_snp_report_get_snps_to_exclude(snp_data1,autosomal_chrs,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep")
excl2 = analyze_snp_report_get_snps_to_exclude(snp_data2,autosomal_chrs,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep")
excl3 = read.table(bad_snps_file,stringsAsFactors = F)[,1]

# Final snp set for subsequent analyses
snp_set = intersect(rownames(snp_data1),rownames(snp_data2))
snp_set = setdiff(snp_set,excl1)
snp_set = setdiff(snp_set,excl2)
snp_set = setdiff(snp_set,excl3)
write.table(t(t(snp_set)),paste(job_dir,"snp_set.txt",sep=""), row.names = F,col.names = F,quote = F)

# take care of the FIDs
fam1 = read.table(paste(input_bfile1,".fam",sep=""),stringsAsFactors = F)
update_ids1 = cbind(fam1[,1:2],paste("file1_",fam1[,1],sep=""),fam1[,2])
write.table(file=paste(job_dir,"update_ids1.txt",sep=""),update_ids1,sep="\t",row.names = F,col.names = F,quote = F)
fam2 = read.table(paste(input_bfile2,".fam",sep=""),stringsAsFactors = F)
# for file 1 we want to remove samples that appear in file 2
inds = !is.element(fam1[,2],set=fam2[,2])
inds = inds & is.element(fam1[,2],set=rownames(sample_metadata_raw))
file1_samples_to_keep = fam1[inds,]
write.table(file=paste(job_dir,"file1_samples_to_keep.txt",sep=""),file1_samples_to_keep,sep="\t",row.names = F,col.names = F,quote = F)
write.table(file=paste(job_dir,"update_ids1.txt",sep=""),update_ids1[inds,],sep="\t",row.names = F,col.names = F,quote = F)

inds2 = is.element(fam2[,2],set=rownames(sample_metadata_raw))
file2_samples_to_keep = fam2[inds2,]
write.table(file=paste(job_dir,"file2_samples_to_keep.txt",sep=""),file2_samples_to_keep,sep="\t",row.names = F,col.names = F,quote = F)

# Correct File 1 before merge
# Remove samples that appear in file 2 (the consortium file)
err_path = paste(job_dir,"prepare_bfile1.err",sep="")
log_path = paste(job_dir,"prepare_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile1,
                 "--keep",paste(job_dir,"file1_samples_to_keep.txt",sep=""),
                 "--extract",paste(job_dir,"snp_set.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "prepare_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

err_path = paste(job_dir,"update_ids_bfile1.err",sep="")
log_path = paste(job_dir,"update_ids_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--update-ids",paste(job_dir,"update_ids1.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "update_ids_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# Correct File 2 before merge
# Remove samples that are not in the metadata
err_path = paste(job_dir,"prepare_bfile2.err",sep="")
log_path = paste(job_dir,"prepare_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile2,
                 "--keep",paste(job_dir,"file2_samples_to_keep.txt",sep=""),
                 "--extract",paste(job_dir,"snp_set.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "prepare_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# Subject QC:
# File 1
err_path = paste(job_dir,"subj_qc_bfile1.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--missing --het",
                 "--out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "subj_qc_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# File 2
err_path = paste(job_dir,"subj_qc_bfile2.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
                 "--missing --het",
                 "--out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "subj_qc_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# Analyze the results
# Subject call rates:
miss1 = read.table(paste(job_dir,"bfile1.imiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.imiss",sep=""),stringsAsFactors = F,header=T)
crs1 = 1-miss1$F_MISS;crs2 = 1-miss2$F_MISS
quantile(crs1)
quantile(crs2)
table(crs1<initial_subj_min_call_rate)
table(crs2<initial_subj_min_call_rate)
het1 = read.table(paste(job_dir,"bfile1.het",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.het",sep=""),stringsAsFactors = F,header=T)
quantile(het1$F)
quantile(het2$F)
table(het2$F < -1,crs2<initial_subj_min_call_rate)
# September 2018: we observed outlier subjects in file 2 only
to_rem = miss2[crs2<initial_subj_min_call_rate,1:2]
remove_subjects_using_plink(paste(job_dir,"bfile2",sep=""),to_rem,job_dir,"file2_initial_subj_qc","bfile2",
                                      batch_script_func=get_sh_default_prefix)

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
wait_for_job()

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
wait_for_job()

# Analyze the SNP QC results
miss1 = read.table(paste(job_dir,"bfile1.lmiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.lmiss",sep=""),stringsAsFactors = F,header=T)

frq1 = read.table(paste(job_dir,"bfile1.frq",sep=""),stringsAsFactors = F,header=T)
frq2 = read.table(paste(job_dir,"bfile2.frq",sep=""),stringsAsFactors = F,header=T)

het1 = read.table(paste(job_dir,"bfile1.hwe",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.hwe",sep=""),stringsAsFactors = F,header=T)

to_rem_file1 = union(miss1$SNP[1-miss1$F_MISS<snp_min_call_rate],
                     frq1$SNP[frq1$MAF < min_maf])
to_rem_file1 = union(to_rem_file1,het1$SNP[het1$P < snp_het_p])

to_rem_file2 = union(miss2$SNP[1-miss2$F_MISS<snp_min_call_rate],
                     frq2$SNP[frq2$MAF < min_maf])
to_rem_file2 = union(to_rem_file2,het2$SNP[het2$P < snp_het_p])

to_rem_files_union = union(to_rem_file1,to_rem_file2)
bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
snps_to_keep = setdiff(bim1[,2],to_rem_files_union)

# Reduce the datasets before the merge
extract_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile1",
    batch_script_func=get_sh_default_prefix)
extract_snps_using_plink(paste(job_dir,"bfile2",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile2",
                         batch_script_func=get_sh_default_prefix)
wait_for_job()

# Read the bim files: we may need to flip some snps
bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
bim2 = read.table(paste(job_dir,"bfile2.bim",sep=""),stringsAsFactors = F)
rownames(bim1) = bim1[,2];rownames(bim2) = bim2[,2]
bim1 = bim1[rownames(bim2),]

get_num_alleles<-function(x){
  return(length(setdiff(unique(x),"0")))
}
xx = cbind(bim1[,5:6],bim2[,5:6])
all_num_alleles = apply(xx,1,get_num_alleles)
one_allele_appears_as_two<-function(x){
  return(get_num_alleles(x)==2 & sum(x=="0")==2)
}
is_one_allele_appears_as_two = apply(xx,1,one_allele_appears_as_two)
snps_to_flip = rownames(bim1)[is_one_allele_appears_as_two | all_num_alleles>2]
length(snps_to_flip)
bim1[snps_to_flip,][1:10,]

flip_nuc<-function(x){
  if(x=="T"){return("A")}
  if(x=="A"){return("T")}
  if(x=="G"){return("C")}
  if(x=="C"){return("G")}
  if(x=="0"){return("0")}
}
flip_snp_info<-function(x){
  return(sapply(x,flip_nuc))
}
repl = apply(xx[snps_to_flip,1:2],1,flip_snp_info)
xx[snps_to_flip,1] = repl[1,]
xx[snps_to_flip,2] = repl[2,]
all_num_alleles2 = apply(xx,1,get_num_alleles)
table(all_num_alleles2)
bim1[all_num_alleles2==0,][1:10,] # should be all NAs because there are no such snps

flip_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_flip,job_dir,"final_snps_to_flip","bfile1",
                         batch_script_func=get_sh_default_prefix)

# Merge
err_path = paste(job_dir,"merge_plink.err",sep="")
log_path = paste(job_dir,"merge_plink.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--bmerge",paste(job_dir,"bfile2",sep=''),
                 "--make-bed --out",paste(job_dir,"merged_mega_data",sep=''))
curr_sh_file = "merge_plink.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

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

flipscan_res = readLines(paste(job_dir,"merged_mega_data.flipscan",sep=""))
arrs = strsplit(flipscan_res[-1],split="\\s+")
names(arrs) = sapply(arrs,function(x)x[3])
table(sapply(arrs,length))
flipscan_failures = sapply(arrs,length) > 11
flipscan_failures = sapply(arrs[flipscan_failures],function(x)x[3])
sapply(arrs[flipscan_failures],function(x)x[2])
bim = read.table(paste(job_dir,"merged_mega_data.bim",sep=""),stringsAsFactors = F)
snps_to_keep = setdiff(bim[,2],flipscan_failures)

extract_snps_using_plink(paste(job_dir,"merged_mega_data",sep=""),snps_to_keep,job_dir,
                         "_final_snps_to_keep_after_flipscan","merged_mega_data",
                         batch_script_func=get_sh_default_prefix)

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # For analysis of specific cohorts
# if(!is.null(analysis_cohorts)){
#   sample_metadata_raw$Cohort = tolower(sample_metadata_raw$Cohort)
#   analysis_cohorts = tolower(analysis_cohorts)
#   curr_samples = rownames(sample_metadata_raw)[is.element(sample_metadata_raw$Cohort,
#                                                           set=analysis_cohorts)]
#   
#   fam_samples = read.table(paste(job_dir,"merged_mega_data.fam",sep=""),stringsAsFactors = F)
#   iid2fid = fam_samples[,1]
#   names(iid2fid) = fam_samples[,2]
#   curr_samples = intersect(curr_samples,fam_samples[,2])
#   curr_samples = cbind(iid2fid[curr_samples],curr_samples)
#   write.table(curr_samples,file=paste(job_dir,"curr_samples.txt",sep=""),
#               row.names = F,quote = F,col.names = F)
#   
#   err_path = paste(job_dir,"keep_cohort_subjects.err",sep="")
#   log_path = paste(job_dir,"keep_cohort_subjects.log",sep="")
#   curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data",sep=''),
#                    "--keep",paste(job_dir,"curr_samples.txt",sep=''),
#                    "--make-bed --out",paste(job_dir,"merged_mega_data",sep=''))
#   curr_sh_file = "keep_cohort_subjects.sh"
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_default_prefix(err_path,log_path),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
#   wait_for_job()
#   system(paste("rm ",job_dir,"merged_mega_data.bed~",sep=""))
#   system(paste("rm ",job_dir,"merged_mega_data.bim~",sep=""))
#   system(paste("rm ",job_dir,"merged_mega_data.fam~",sep=""))
# }
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# impute sex: eahc of the Mega datasets separately
run_plink_sex_check_x_chrom<-function(bfile,out_path,
                                      batch_script_func=get_sh_default_prefix,...){
  err_path = paste(job_dir,"impute_sex.err",sep="")
  log_path = paste(job_dir,"impute_sex.log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--check-sex --chr 22-24 --out",bfile)
  curr_sh_file = paste(out_path,"impute_sex.sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
run_plink_sex_check_x_chrom(input_bfile1,job_dir)
run_plink_sex_check_x_chrom(input_bfile2,job_dir)

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
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 chr_filter,
                 "--make-bed --out",paste(job_dir,"merged_mega_data_autosomal",sep=''))
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
# # Exclude low maf snps.
# # Also add freq, pca, and misingness analyses
# if(is.null(analysis_cohorts)){
#   maf_snps_to_exclude = snp_data$Name[snp_data$Minor.Freq<min_maf]
#   write.table(t(t(as.character(snp_data$Name[maf_snps_to_exclude]))),
#               file=paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
#               row.names = F,col.names = F,quote = F)
#   err_path = paste(job_dir,"maf_filter.err",sep="")
#   log_path = paste(job_dir,"maf_filter.log",sep="")
#   curr_cmd = paste("plink --bfile",paste(job_dir,"chr_filter",sep=''),
#                    "--exclude",paste(job_dir,"maf_snps_to_exclude.txt",sep=''),
#                    "--maf",min_maf,
#                    "--pca --freq --missing",
#                    "--make-bed --out",paste(job_dir,"maf_filter_data",sep=''))
#   curr_sh_file = "maf_filter.sh"
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_default_prefix(err_path,log_path),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
#   wait_for_job()
# }

# if(!is.null(analysis_cohorts)){
#   err_path = paste(job_dir,"maf_filter.err",sep="")
#   log_path = paste(job_dir,"maf_filter.log",sep="")
#   curr_cmd = paste("plink --bfile",paste(job_dir,"chr_filter",sep=''),
#                    "--maf 0.01", # for elite, consider changing for other cohorts
#                    "--pca --freq --missing",
#                    "--make-bed --out",paste(job_dir,"maf_filter_data",sep=''))
#   curr_sh_file = "maf_filter.sh"
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_default_prefix(err_path,log_path),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
#   wait_for_job()
# }

####################################################################################################
####################################################################################################
####################################################################################################
# Here we analyze the results of some of the jobs above

# compare to known sex from the metadata
metadata_sex = sample_metadata_raw$Sex_input_data
names(metadata_sex) = rownames(sample_metadata_raw)

imputed_sex1 = read_plink_table(paste(input_bfile1,".sexcheck",sep=''))[,4]
imputed_sex2 = read_plink_table(paste(input_bfile2,".sexcheck",sep=''))[,4]

inds1 = intersect(names(metadata_sex),names(imputed_sex1))
table(metadata_sex[inds1],imputed_sex1[inds1])
inds2 = intersect(names(metadata_sex),names(imputed_sex2))
table(metadata_sex[inds2],imputed_sex2[inds2])

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

