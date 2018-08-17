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

####################################################################################################
####################################################################################################
####################################################################################################
# Prepare the GWAS pheno data

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

# write phe files
pheno_data = cbind(as.character(iid_to_fid[remaining_samples]),pheno_data)
colnames(pheno_data) = c("FID","IID","sex","ExerciseGroup","Batch","Age",paste("PC",1:10,sep=""))
table(pheno_data$Age)
write.table(file=paste(job_dir,"three_group_analysis_genepool_controls.phe",sep=''),
            pheno_data,sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=paste(job_dir,"three_group_analysis_sex_update.phe",sep=''),
            pheno_data[,1:3],sep=" ",row.names = F,col.names = T,quote=F)


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

# 3. Logistic Cooper vs. Genepool without age
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

# 4. Logistic Cooper vs. Elite with age
table(pheno_data$ExerciseGroup)
sample_inds = pheno_data$ExerciseGroup != "-1"
curr_pheno = pheno_data[sample_inds,]
curr_pheno[,4] = as.character(as.numeric(curr_pheno[,4]+1))
pheno_file = paste(job_dir,"cooper_vs_elite.phe",sep='')
write.table(file=pheno_file,curr_pheno[,c(1:2,4)],sep=" ",row.names = F,col.names = T,quote=F)
covar_file = paste(job_dir,"cooper_vs_elite.phe",sep='')
write.table(file=covar_file,curr_pheno[,-4],sep=" ",row.names = F,col.names = T,quote=F)
#jobs_before = get_my_jobs()
err_path = paste(job_dir,"cooper_vs_elite.err",sep="")
log_path = paste(job_dir,"cooper_vs_elite.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--logistic hide-covar firth-fallback",
                 "--1",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Age,Batch,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"cooper_vs_elite",sep=''))
curr_sh_file = "cooper_vs_elite.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Repeat the GWAS above after applying PCA and clustering, number of clusters is a parameter
# given by the user. This is determined based on either some prior knowledge or by looking at the
# elbow plot of kmeans.
d = covariate_matrix
d2 = sample_metadata_raw
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids

# Cluster by the first two PCs
pc_x = as.matrix(d[,c("PC1","PC2")])
pc_x_kmeans = kmeans(pc_x,num_pca_clusters)
table(pc_x_kmeans$cluster)
table(pc_x_kmeans$cluster,d$Cohort)
kmeans_res = pc_x_kmeans$cluster

# Select subjects from the largest cluster
clustable = table(kmeans_res)
selected_cluster = names(which(clustable == max(clustable)))
selected_subjects = names(which(kmeans_res == selected_cluster))

# 0. Linear for all groups + age, sex, 5 PCs
sample_inds = is.element(rownames(pheno_data),set=selected_subjects)
pheno_file = paste(job_dir,"three_group_analysis_genepool_controls_pcafilter.phe",sep='')
covar_file = paste(job_dir,"three_group_analysis_genepool_controls_pcafilter.phe",sep='')
curr_pheno = pheno_data[sample_inds,]
write.table(file=pheno_file,curr_pheno[,c(1:2,4)],sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=covar_file,curr_pheno[,-4],sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(job_dir,"gwas_three_groups_linear_pcafilter.err",sep="")
log_path = paste(job_dir,"gwas_three_groups_linear_pcafilter.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"gwas_three_groups_linear_pcafilter",sep=''))
curr_sh_file = "gwas_three_groups_linear_pcafilter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# 4. Logistic Cooper vs. Elite with age
sample_inds = is.element(rownames(pheno_data),set=selected_subjects) & pheno_data$ExerciseGroup != "-1"
curr_pheno = pheno_data[sample_inds,]
pheno_file = paste(job_dir,"cooper_vs_elite_pcafilter.phe",sep='')
write.table(file=pheno_file,curr_pheno[,c(1:2,4)],sep=" ",row.names = F,col.names = T,quote=F)
covar_file = paste(job_dir,"cooper_vs_elite_covars_pcafilter.phe",sep='')
write.table(file=covar_file,curr_pheno[,-4],sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(job_dir,"cooper_vs_elite_pcafilter.err",sep="")
log_path = paste(job_dir,"cooper_vs_elite_pcafilter.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"maf_filter",sep=''),
                 "--logistic hide-covar firth-fallback",
                 "--1",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Age,Batch,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(job_dir,"cooper_vs_elite_pcafilter",sep=''))
curr_sh_file = "cooper_vs_elite_pcafilter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Print all GWAS results in FUMA's format
create_fuma_files_for_fir(job_dir,
                          paste(job_dir,"maf_filter.bim",sep=""),
                          paste(job_dir,"maf_filter.frq",sep=""))








