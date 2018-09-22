# This script handles the GWAS analysis of our dataset 
# i.e., all analyses that do not require external controls
# Currently: not sure this is necessary given that we omit genepool for now

####################################################################################################
####################################################################################################
####################################################################################################

# Define Input parameters
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/"
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

####################################################################################################
####################################################################################################
####################################################################################################
# A set of GWAS runs starts here

# # 0. Linear for all groups + age, sex, 5 PCs
# pheno_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
# covar_file = paste(job_dir,"three_group_analysis_genepool_controls.phe",sep='')
# #jobs_before = get_my_jobs()
# err_path = paste(job_dir,"gwas_three_groups_linear.err",sep="")
# log_path = paste(job_dir,"gwas_three_groups_linear.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,"maf_filter",sep=''),
#                  "--linear hide-covar",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name ExerciseGroup"),
#                  "--allow-no-sex",
#                  paste("--covar",covar_file),
#                  "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
#                  "--adjust",
#                  "--out",paste(job_dir,"gwas_three_groups_linear",sep=''))
# curr_sh_file = "gwas_three_groups_linear.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()

# # 3. Logistic Cooper vs. Genepool without age
# table(pheno_data$ExerciseGroup)
# sample_inds = pheno_data$ExerciseGroup != "1"
# curr_pheno = pheno_data[sample_inds,]
# curr_pheno[,4] = as.character(as.numeric(curr_pheno[,4]+1))
# pheno_file = paste(job_dir,"cooper_vs_genepool.phe",sep='')
# write.table(file=pheno_file,curr_pheno[,c(1:2,4)],sep=" ",row.names = F,col.names = T,quote=F)
# covar_file = paste(job_dir,"cooper_vs_genepool_covar.phe",sep='')
# write.table(file=covar_file,curr_pheno[,-4],sep=" ",row.names = F,col.names = T,quote=F)
# #jobs_before = get_my_jobs()
# err_path = paste(job_dir,"cooper_vs_genepool.err",sep="")
# log_path = paste(job_dir,"cooper_vs_genepool.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,"maf_filter",sep=''),
#                  "--logistic hide-covar firth-fallback",
#                  "--1",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name ExerciseGroup"),
#                  "--allow-no-sex",
#                  paste("--covar",covar_file),
#                  "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5",
#                  "--adjust",
#                  "--out",paste(job_dir,"cooper_vs_genepool",sep=''))
# curr_sh_file = "cooper_vs_genepool.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()

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
pc_x_kmeans = kmeans(pc_x,5)
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








