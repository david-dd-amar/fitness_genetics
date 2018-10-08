
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# assumption: merged bed has frq and pca results
# # all cohorts together
# bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/merged_data_qctool_bed"
# external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/new_bed_1.frq"
# our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/new_bed_2.frq"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/gwas/"
# 
# # elite and ukbb alone
# bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/with_ukbb/merged_data_qctool_bed"
# external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/with_ukbb/new_bed_1.frq"
# our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/with_ukbb/new_bed_2.frq"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/with_ukbb/gwas/"

# September 2018: new MEGA analysis, HRC as the panel
bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/merged_data_qctool_bed"
alt_bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/merged_data_plink"
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/new_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/new_bed_2.frq"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/gwas/"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/integrated_sample_metadata_and_covariates.phe"

# September 2018: new MEGA analysis, 1000g as the panel
bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/merged_data_qctool_bed"
alt_bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/merged_data_plink"
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/new_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/new_bed_2.frq"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/gwas/"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/integrated_sample_metadata_and_covariates.phe"

# September 2018: new MEGA analysis, 1000g as the panel (second run, better preprocessing)
bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/merged_data_qctool_bed"
alt_bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/merged_data_plink"
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/new_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/new_bed_2.frq"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/gwas/"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/integrated_sample_metadata_and_covariates_after_pca1.phe"

external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"
our_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"

try({system(paste("mkdir",out_path))})

our_data_mafs_by_group = list(
  "genepool" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/genepool_cohort_freq.frq",
  "cooper" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/Cooper_cohort_freq.frq",
  "elite" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/ELITE_cohort_freq.frq"
)

# The steps of the analysis below
# 1. Use our dataset to determine "white" subjects
# 2. Make sure our "white" subjects are well clustered in the new pca
# 3. If (2) is successfull: select N controls that are close to our ELITE and Cooper samples
# 4. Also define the white genepool controls
# 5. Run GWASs
# We also do some QC analyses locally below. 

####################################################################################################
####################################################################################################
####################################################################################################
# Read covars, pca, and create phe file
# Our covars file has the exercise group in column 4
# This phenotype is encoded as 1 = elite, 0=cooper, and -1=genepool
# Other covariates in this file are: sex, age, batch, pcs
# First two columns are FID and IID
# This file is one of the results of the gwas_flow.R script, the gwas run
# part of this file should be revised (June 2018)
our_covars = read.table(our_covars_path,header=T,stringsAsFactors = F)
our_phe = as.character(our_covars[,"Cohort"])
table(our_covars[,"Cohort"])
cohorts = as.character(our_covars[,"Cohort"])
our_covars[,"Cohort"] = our_covars[,"Cohort"] + 1

# Read external DB info
external_covars = read.table(external_covars_path,stringsAsFactors = F)
external_covars = cbind(as.character(external_covars[,1]),external_covars)
external_samples = as.character(external_covars[,1])

# For ukbb - add batches
batch_data = as.matrix(
  read.table("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.tab.fam",
             stringsAsFactors = F, header=T))
rownames(batch_data) = as.character(batch_data[,1])
external_covars = cbind(batch_data[external_samples,c(1:2,5:6)],external_covars[,4])
external_covars = cbind(external_covars,rep("1",nrow(external_covars)),rep("ukbb",nrow(external_covars)))
colnames(external_covars) = c("FID","IID","sex","Batch","Age","ExerciseGroup","CohortName")

# Define the final covariance matrix (with the ExerciseGroup column), without the PCs
our_covars_wo_pcs = cbind(our_covars[,c("FID","IID","sex","batch","age","Cohort")],cohorts)
colnames(our_covars_wo_pcs) = colnames(external_covars)
covars = as.matrix(rbind(our_covars_wo_pcs,external_covars))
rownames(covars) = covars[,"IID"]
covars[covars[,7]=="2",7] = "elite"
covars[covars[,7]=="1",7] = "cooper"

# Define the PCAs: the combined dataset and ours
combined_pcs = read_pca_res(paste(bfile,".eigenvec",sep=""))
our_dataset_pcs = our_covars[,grepl("^PC",colnames(our_covars))]
rownames(our_dataset_pcs) = our_covars[,"IID"]
pcs_explained_var = read.table(paste(bfile,".eigenval",sep=""))

# Compare the pcs of the bed and the alternative bed
combined_pcs = read_pca_res(paste(bfile,".eigenvec",sep=""))
combined_pcs2 = read_pca_res(paste(alt_bfile,".eigenvec",sep=""))
length(intersect(rownames(combined_pcs),rownames(combined_pcs2))) == nrow(combined_pcs)
combined_pcs2 = combined_pcs2[rownames(combined_pcs),1:40]
pc_corrs = cor(combined_pcs,combined_pcs2)
apply(abs(pc_corrs),1,max)
apply(abs(pc_corrs),2,max)
pc_corrs2 = cor(combined_pcs[rownames(our_dataset_pcs),],our_dataset_pcs)
apply(abs(pc_corrs2),1,max)
apply(abs(pc_corrs2),2,max)

subjects_for_analysis = intersect(covars[,"IID"],rownames(combined_pcs))
covars = cbind(covars[subjects_for_analysis,],combined_pcs[subjects_for_analysis,])
# covars[,"Batch"] = cov_phe_col_to_plink_numeric_format(covars[,"Batch"])
for(j in 1:ncol(covars)){
  covars[,j] = gsub(" ","",as.character(covars[,j]))
}

write.table(file=paste(out_path,"all_cohorts.phe",sep=''),
            covars,sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# Preprocessing the subject set: relatedness and clustering

# Analyze the relatedness report
library("igraph",lib.loc = "~/R/packages")
rl_data = read.table(paste(bfile,".genome",sep=""),header=T,stringsAsFactors = F)
rl_edges = as.matrix(rl_data[,c("IID1","IID2")])
mode(rl_edges) = "character"
rl_g = igraph::graph_from_edgelist(rl_edges,directed = F)
rl_clusters = clusters(rl_g)[[1]]
rl_subjects_to_remove = c()
for(cl in unique(rl_clusters)){
  curr_subjects = names(rl_clusters)[rl_clusters==cl]
  rl_subjects_to_remove = c(rl_subjects_to_remove,curr_subjects[-1])
}

# Cluster the data, take the largest cluster and continue for the analysis
d = read.table(paste(out_path,"all_cohorts.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
d2 = read.delim(our_metadata,stringsAsFactors = F)
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID;altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids; names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids; table(is_jap)
cohorts = d$CohortName; table(cohorts)
d2_analysis_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
jap_samples = d2_analysis_ids [is_jap]
alldata_is_jap = is.element(d$IID,set=jap_samples)
names(alldata_is_jap) = d$IID

set.seed(123)
pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
rownames(pc_x) = rownames(d)
dd = dist(pc_x,method="manhattan")
h = hclust(dd,method = "single")
kmeans_res = run_hclust(pc_x,150,dd,h)
kmeans_res[kmeans_res!=1] = 0
table(kmeans_res)
table(kmeans_res,d$Cohort)

to_rem = rep(F,nrow(d))
for(j in 1:20){
  x = d[,paste("PC",j,sep="")]
  x = (x-mean(x))/sd(x)
  print(sum(abs(x)>8))
  to_rem[abs(x)>8] = T
}
table(to_rem,d$CohortName)
table(to_rem,kmeans_res)
table(kmeans_res[!to_rem],d$CohortName[!to_rem])

kmeans_res[to_rem] = 0
kmeans_res[rl_subjects_to_remove] = 0
table(kmeans_res)
table(kmeans_res,d$Cohort)

# set.seed(123)
# pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
# rownames(pc_x) = d$IID
# wss <- sapply(1:10,
#               function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
# wss[2:length(wss)]/wss[1:(length(wss)-1)]
# ## Kmeans-based analysis
# kmeans_res <- kmeans(pc_x, 5)$cluster
# table(kmeans_res)
# table(kmeans_res,d[rownames(pc_x),]$CohortName)
# table(kmeans_res,alldata_is_jap[rownames(pc_x)]) # Japanese are well clustered and removed

pc_ps = c()
for(j in 1:40){
  pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
    pc = d[,paste("PC",j,sep="")],
    y = d[,"CohortName"]
  )
}
pc_ps = p.adjust(pc_ps)
pc_inds = which(pc_ps < 0.01) # Before correction: almost all

# get the largest cluster and take its subjects
cl_tb = table(kmeans_res)
cl_fa = names(cl_tb)[cl_tb==max(cl_tb)]
selected_subjects_for_gwas = names(kmeans_res)[kmeans_res==cl_fa]
selected_subjects_for_gwas = setdiff(selected_subjects_for_gwas,rl_subjects_to_remove)
d = d[selected_subjects_for_gwas,]

dim(d)
save(selected_subjects_for_gwas,kmeans_res,pc_x,file=paste(out_path,"clustering_data.RData",sep=""))
write.table(file=paste(out_path,"kmeans_cleaned.phe",sep=''),
            d,sep=" ",row.names = F,col.names = T,quote=F)

print(paste("After clustering and relatedness analysis, number of remaining samples:",
            length(selected_subjects_for_gwas)))
curr_fam = read.table(paste(bfile,".fam",sep=""),stringsAsFactors = F,header = F)
rownames(curr_fam) = as.character(curr_fam[,2])
curr_fam = curr_fam[setdiff(rownames(curr_fam),selected_subjects_for_gwas),1:2]
remove_subjects_using_plink(bfile,curr_fam,
                            out_path,"_pca_and_rl_subj_qc","merged_data_after_pca_rl_filters",
                            batch_script_func=get_sh_default_prefix)
wait_for_job()
print("After PCA and Rl analysis, data sizes are:")
print(paste("number of samples:",length(readLines(paste(out_path,"merged_data_after_pca_rl_filters.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(out_path,"merged_data_after_pca_rl_filters.bim",sep="")))))

####################################################################################################
####################################################################################################
####################################################################################################

# Snp prune
analysis_name = "plink_prune"
err_path = paste(out_path,analysis_name,"_ld_report.err",sep="")
log_path = paste(out_path,analysis_name,"_ld_report.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_after_pca_rl_filters",sep=''),
                 "--indep-pairwise 250 10",0.1,
                 "--maf 0.01",
                 "--out",paste(out_path,analysis_name,sep=""))
curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()
# Run PCA
err_path = paste(out_path,"run_pca.err",sep="")
log_path = paste(out_path,"sun_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_after_pca_rl_filters",sep=''),
                 "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
                 "--pca 40 --out",paste(out_path,"merged_data_after_pca_rl_filters",sep=''))
curr_sh_file = "run_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# Freq on all snps
err_path = paste(out_path,"run_frq.err",sep="")
log_path = paste(out_path,"run_frq.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_after_pca_rl_filters",sep=''),
                 "--freq --out",paste(out_path,"merged_data_after_pca_rl_filters",sep=''))
curr_sh_file = "run_frq.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

# update our covariates
new_pca_res = read_pca_res(paste(out_path,"merged_data_after_pca_rl_filters.eigenvec",sep=""))
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
d = d[rownames(new_pca_res),]
d[,colnames(new_pca_res)] = new_pca_res
write.table(d,file=paste(out_path,"kmeans_cleaned.txt",sep=''),
            sep="\t",quote=F,row.names = F)
write.table(d,file=paste(out_path,"kmeans_cleaned.phe",sep=''),
            sep=" ",quote=F,row.names = F)

# Compute the new PCs association
pc_ps = c()
for(j in 1:40){
  pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
    pc = d[,paste("PC",j,sep="")],
    y = d[,"CohortName"]
  )
}
pc_ps = p.adjust(pc_ps)
pc_inds = which(pc_ps < 0.001) # Before correction: almost all
PCs = paste("PC",pc_inds[pc_inds<20],sep="")

# # Take the closest sample to each of our subjects and recalculate
# our_samples = rownames(d)[d$CohortName != "ukbb"]
# ukbb_samples = rownames(d)[d$CohortName == "ukbb"]
# pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
# rownames(pc_x) = rownames(d)
# selected_samples = c()
# n_to_select = 2
# for (i in our_samples){
#   curr_dists = sweep(pc_x[ukbb_samples,],2,pc_x[i,])
#   curr_dists = curr_dists^2
#   curr_dists = sqrt(rowSums(curr_dists))
#   curr_dists = sort(curr_dists,decreasing = F)
#   selected_samples = union(selected_samples,names(curr_dists)[1:n_to_select])
#   print(length(selected_samples))
# }
# subjects_for_analysis = c(our_samples,selected_samples)
# pc_ps = c()
# for(j in 1:40){
#   pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
#     pc = d[subjects_for_analysis,paste("PC",j,sep="")],
#     y = d[subjects_for_analysis,"CohortName"]
#   )
# }
# pc_ps = p.adjust(pc_ps)
# pc_inds = which(pc_ps < 0.01) # Before correction: almost all

####################################################################################################
####################################################################################################
####################################################################################################

# Run GWAS
covar_file = paste(out_path,"kmeans_cleaned.phe",sep='')
gwas_bfile = paste(out_path,"merged_data_after_pca_rl_filters",sep='')

# 1. Linear of all three groups + sex, age, and up to 20 PCs
err_path = paste(out_path,"gwas_three_groups_linear.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",gwas_bfile,"--linear hide-covar",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 paste("--covar-name sex,Age,",paste(PCs,collapse=","),sep=""),
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 2. Logistic: Elite vs. UKBB, + sex, age, and up to 20 PCs
covars_copy = d[d$CohortName!="cooper",]
covars_copy$ExerciseGroup[covars_copy$ExerciseGroup=="3"] = 2
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_elite.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_elite.err",sep="")
log_path = paste(out_path,"ukbb_vs_elite.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 paste("--covar-name sex,Age,",paste(PCs,collapse=","),sep=""),
                 "--adjust --out",paste(out_path,"ukbb_vs_elite_logistic",sep=''))
curr_sh_file = "ukbb_vs_elite_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 3. Logistic: Cooper vs. UKBB, + sex, age, and up to 20 PCs
covars_copy = d[d$CohortName!="elite",]
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_cooper.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_cooper.err",sep="")
log_path = paste(out_path,"ukbb_vs_cooper.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 paste("--covar-name sex,Age,",paste(PCs,collapse=","),sep=""),
                 "--adjust --out",paste(out_path,"ukbb_vs_cooper_logistic",sep=''))
curr_sh_file = "ukbb_vs_cooper_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################

d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID

# Take the closest sample to each of our subjects and recalculate
our_samples = rownames(d)[d$CohortName == "elite"]
ukbb_samples = rownames(d)[d$CohortName == "ukbb"]
pc_x = as.matrix(d[,paste("PC",1:10,sep="")])
rownames(pc_x) = rownames(d)
selected_samples = c()
n_to_select = 2
for (i in our_samples){
  curr_dists = sweep(pc_x[ukbb_samples,],2,pc_x[i,])
  curr_dists = curr_dists^2
  curr_dists = sqrt(rowSums(curr_dists))
  curr_dists = sort(curr_dists,decreasing = F)
  selected_samples = union(selected_samples,names(curr_dists)[1:n_to_select])
  print(length(selected_samples))
}
subjects_for_analysis = c(our_samples,selected_samples)
pc_ps = c()
for(j in 1:40){
  pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
    pc = d[subjects_for_analysis,paste("PC",j,sep="")],
    y = d[subjects_for_analysis,"CohortName"]
  )
}
pc_ps = p.adjust(pc_ps)
pc_inds = which(pc_ps < 0.01) # Before correction: almost all

# 4. Logistic: Elite vs. UKBB, + sex, age, and 10 PCs
covars_copy = d[subjects_for_analysis,]
covars_copy$ExerciseGroup[covars_copy$ExerciseGroup=="3"] = 2
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_elite_matched_samples.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_elite_matched_samples.err",sep="")
log_path = paste(out_path,"ukbb_vs_elite_matched_samples.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 paste("--covar-name sex,Age,",paste(paste("PC",1:10,sep=""),collapse=","),sep=""),
                 "--adjust --out",paste(out_path,"ukbb_vs_elite_logistic_matched_samples",sep=''))
curr_sh_file = "ukbb_vs_elite_logistic_matched_samples.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# # 5. Logistic: Cooper vs. UKBB, + sex, age, and 15 PCs
# covars_copy = d[d$CohortName!="elite",]
# table(covars_copy$ExerciseGroup)
# covar_file = paste(out_path,"ukbb_vs_cooper.phe",sep='')
# write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
# err_path = paste(out_path,"ukbb_vs_cooper.err",sep="")
# log_path = paste(out_path,"ukbb_vs_cooper.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
#                  paste("--pheno",covar_file),
#                  paste("--pheno-name ExerciseGroup"),
#                  paste("--covar",covar_file),
#                  paste("--covar-name sex,Age,",paste(paste("PC",1:15,sep=""),collapse=","),sep=""),
#                  "--adjust --out",paste(out_path,"ukbb_vs_cooper_logistic",sep=''))
# curr_sh_file = "ukbb_vs_cooper_logistic.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# # Sanity check 1 : UKBB vs. Genepool
# load(paste(out_path,"clustering_data.RData",sep=""))
# d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
# rownames(d) = d$IID
# d = d[selected_subjects_for_gwas,]
# covars_copy = d[d$CohortName!="elite" & d$CohortName!="cooper",]
# covars_copy$ExerciseGroup[covars_copy$CohortName=="ukbb"]="1"
# covars_copy$ExerciseGroup[covars_copy$CohortName=="genepool"]="2"
# table(covars_copy$ExerciseGroup)
# covar_file = paste(out_path,"ukbb_vs_genepool.phe",sep='')
# write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
# err_path = paste(out_path,"ukbb_vs_genepool.err",sep="")
# log_path = paste(out_path,"ukbb_vs_genepool.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
#                  paste("--pheno",covar_file),
#                  paste("--pheno-name ExerciseGroup"),
#                  paste("--covar",covar_file),
#                  "--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
#                  "--adjust --out",paste(out_path,"ukbb_vs_genepool_logistic",sep=''))
# curr_sh_file = "ukbb_vs_genepool_logistic.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Sanity check 2: UKBB vs. not UKBB: flip scan
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
covars_copy = d
covars_copy$ExerciseGroup[covars_copy$CohortName=="ukbb"]="1"
covars_copy$ExerciseGroup[covars_copy$CohortName!="ukbb"]="2"
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_nonukbb.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_nonukbb.err",sep="")
log_path = paste(out_path,"ukbb_vs_nonukbb.log",sep="")
curr_cmd = paste("plink",
                 "--bfile",gwas_bfile,"--logistic --flip-scan --allow-no-sex --test-missing",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--adjust --out",paste(out_path,"ukbb_vs_nonukbb_logistic",sep=''))
curr_sh_file = "ukbb_vs_nonukbb_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size = 10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################
# Run GWAS for each PC
covar_file = paste(out_path,"kmeans_cleaned.phe",sep='')
for (j in 1:20){
  err_path = paste(out_path,"gwas_PC",j,".err",sep="")
  log_path = paste(out_path,"gwas_PC",j,".log",sep="")
  curr_cmd = paste("plink2",
                   "--bfile",gwas_bfile,"--linear hide-covar",
                   paste("--pheno-name",paste("PC",j,sep="")),
                   paste("--pheno",covar_file),
                   paste("--covar",covar_file),
                   "--covar-name sex,Age",
                   "--adjust",
                   "--out",paste(out_path,"gwas_PC",j,"",sep=''))
  curr_sh_file = paste("gwas_PC",j,".sh",sep="")
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}

####################################################################################################
####################################################################################################
####################################################################################################
# output all gwas results into a new dir with input files for fuma interpretation

create_fuma_files_for_fir(out_path,
                          paste(gwas_bfile,".bim",sep=""),
                          paste(gwas_bfile,".frq",sep=""),p = 1,maf = 0.05,
                          snps_to_exclude_from_results=NULL)

# Compare the results
res_files = list.files(out_path)
res_files = res_files[grepl("adjusted$",res_files)]
m =  NULL
for (f in res_files){
  res = read.table(paste(out_path,f,sep=''),stringsAsFactors = F)
  p = res[,3];names(p) = res[,2]
  if(is.null(m)){m=p;next}
  if(is.null(dim(m))){m = cbind(m,p[names(m)]);next}
  m = cbind(m,p[rownames(m)])
}
table(m[,1]<5e-8,m[,2]<5e-8)
table(m[,3]<0.001,m[,2]<5e-8)

# Optional: discard unreliable UKBB snps
bad_ukbb_snps = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/bad_ukbb_snps.txt",
                           stringsAsFactors = F)
bad_ukbb_snps = bad_ukbb_snps[,1]
create_fuma_files_for_fir(out_path,
                          paste(gwas_bfile,".bim",sep=""),
                          paste(gwas_bfile,".frq",sep=""),p = 1,maf = 0.01,
                          snps_to_exclude_from_results=bad_ukbb_snps)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Analysis of the sanity check results

# # mafs
# get_mafs_from_file<-function(path){
#   freqs = read.table(path,header=T,stringsAsFactors = F)
#   rownames(freqs) = freqs[,2]
#   mafs = freqs$MAF;names(mafs)=rownames(freqs)
#   return(mafs)
# }
# combined_mafs = get_mafs_from_file(paste(bfile,".frq",sep=""))
# our_mafs = get_mafs_from_file(our_data_mafs)
# external_mafs = get_mafs_from_file(external_data_mafs)
# names(our_mafs) = names(external_mafs)
# our_group_mafs = lapply(our_data_mafs_by_group,get_mafs_from_file)
# cor(our_group_mafs[[1]],our_group_mafs[[2]])
# 
# # gwas'
# ukbb_vs_gp_res = read.table(
#   paste(out_path,"ukbb_vs_genepool_logistic.ExerciseGroup.glm.logistic.hybrid.adjusted",sep=''),
#   header=F,stringsAsFactors = F,check.names = F)
# rownames(ukbb_vs_gp_res) = ukbb_vs_gp_res[,2]
# colnames(ukbb_vs_gp_res) = ukbb_vs_gp_res[1,]
# ukbb_vs_gp_res = ukbb_vs_gp_res[-1,]
# ukbb_vs_nonukbb_res = read.table(
#   paste(out_path,"ukbb_vs_nonukbb_logistic.missing.adjusted",sep=""),
#   header=F,stringsAsFactors = F
# )
# rownames(ukbb_vs_nonukbb_res) = ukbb_vs_nonukbb_res[,2]
# colnames(ukbb_vs_nonukbb_res) = ukbb_vs_nonukbb_res[1,]
# ukbb_vs_nonukbb_res = ukbb_vs_nonukbb_res[-1,]

flipscan_res = read.delim2(
  paste(out_path,"ukbb_vs_nonukbb_logistic.flipscan",sep=""),
  header=T,stringsAsFactors = F,na.strings = NULL,sep="\t"
)
# to interpret these results see: http://zzz.bwh.harvard.edu/plink/dataman.shtml#flipscan
# basically: snps whose num negatives (column 9) is > 0 are problematic
flipscan_res = apply(flipscan_res,1,function(x)strsplit(x,split="\\s+")[[1]])
flengths = sapply(flipscan_res,length)
negs = as.numeric(sapply(flipscan_res,function(x)x[10]))
table(flengths==12,negs>0)
flipscan_problematic_snps = sapply(flipscan_res[flengths==12],function(x)x[3])

# Look at the results
ps_check1 = as.numeric(ukbb_vs_gp_res[,ncol(ukbb_vs_gp_res)])
names(ps_check1) = rownames(ukbb_vs_gp_res)
ps_check2 = as.numeric(ukbb_vs_nonukbb_res[,ncol(ukbb_vs_nonukbb_res)])
names(ps_check2) = rownames(ukbb_vs_nonukbb_res)
# inds = intersect(names(ps_check2),names(ps_check1))
# cor(ps_check1[inds],ps_check2[inds])
check1_snps = names(ps_check1)[ps_check1<0.00001]
check2_snps = names(ps_check2)[ps_check2<0.00001] # means that missing values are not random

# topX=10000
# tokeep = rep(T,length(our_mafs))
# names(tokeep)=names(our_mafs)
# while(sum(tokeep)>100000){
#   x = our_mafs[tokeep]
#   y = external_mafs[tokeep]
#   l = lm(y~x)
#   print(summary(l))
#   r = residuals(l)
#   absr = abs(r)
#   currthr = sort(absr,decreasing = T)[topX]
#   toremove = names(absr)[absr>=currthr]
#   tokeep[toremove] = F
# }

# table(mafs1[check1_snps]<0.01,mafs2[check1_snps]<0.01)
# length(intersect(flipscan_problematic_snps,check1_snps))

# snps to ignore:
low_mafs1 = names(our_mafs)[our_mafs < 0.01]
low_mafs2 = names(external_mafs)[external_mafs<0.01]
snps_to_exclude_from_results = union(flipscan_problematic_snps,low_mafs1)
snps_to_exclude_from_results = union(snps_to_exclude_from_results,low_mafs2)

# Are these the zero p-value snps in our initial analysis?
gwas_res_example = read.table(paste(out_path,"gwas_three_groups_linear.ExerciseGroup.glm.linear.adjusted",sep=""),
                              stringsAsFactors = F)
zero_pval_snps = gwas_res_example[gwas_res_example[,ncol(gwas_res_example)] < 1e-50 , 2]
table(is.element(zero_pval_snps,set=snps_to_exclude_from_results))
table(is.element(check1_snps,set=snps_to_exclude_from_results))

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Print all GWAS results in FUMA's format
create_fuma_files_for_fir(out_path,
  paste(bfile,".bim",sep=""),
  paste(bfile,".frq",sep=""),p = 1,maf = 0.01,
  snps_to_exclude_from_results=NULL)

# Test: elite vs. ukbb: linear vs. logistic (same analysis basically)
gwas_res_example1 = read.table(paste(out_path,"gwas_three_groups_linear.ExerciseGroup.glm.linear.adjusted",sep=""),
                              stringsAsFactors = F)
gwas_res_example2 = read.table(paste(out_path,"ukbb_vs_elite_logistic.ExerciseGroup.glm.logistic.hybrid.adjusted",sep=""),
                               stringsAsFactors = F)
rownames(gwas_res_example1) = gwas_res_example1[,2]
rownames(gwas_res_example2) = gwas_res_example2[,2]
setdiff(gwas_res_example1[,2],gwas_res_example2[,2])
gwas_res_example2 = gwas_res_example2[gwas_res_example1[,2],]
ps1 = gwas_res_example1[,10];ps2=gwas_res_example2[,10]
ps1[ps1==0] = 1e-200
ps2[ps2==0] = 1e-200
cor(log(ps1),log(ps2))
table(ps1<1e-50,ps2<1e-8)
low_ps1_snps = ps1<1e-50
low_ps2_snps = ps2<1e-8

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

setwd("/Users/David/Desktop/elite/sept2018_prepro_res/with_ukbb/")
d = read.table("all_cohorts.phe",header=T,stringsAsFactors = F)
# d = read.table("kmeans_cleaned.phe",header=T,stringsAsFactors = F)
rownames(d) = d$IID
d2 = read.delim("../../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt",stringsAsFactors = F)
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID;altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids; names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids; table(is_jap)
cohorts = d$CohortName; table(cohorts)
d2_analysis_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
jap_samples = d2_analysis_ids [is_jap]
alldata_is_jap = is.element(d$IID,set=jap_samples)
names(alldata_is_jap) = d$IID
# our_pca = read_pca_res("../../analysis/final_dataset_for_analysis.eigenvec")

# # Examine the PCA results
# library(corrplot)
# pca1 = read_pca_res("merged_data_plink.eigenvec")
# pca2 = read_pca_res("merged_data_qctool_bed.eigenvec")
# pca2 = pca2[rownames(pca1),]
# all(rownames(pca1)==rownames(pca2))
# corrs = cor(pca1,pca2)
# corrplot(corrs)
# pcainds = intersect(rownames(pca1),rownames(our_pca))
# corrs = cor(pca1[pcainds,],our_pca[pcainds,])
# corrplot(corrs)
# d = d[rownames(pca1),]

# Cluster by PCs
pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
# pcs_explained_var = read.table("merged_data_qctool_bed.eigenval")[,1]
# for(j in 1:ncol(pc_x)){pc_x[,j]=pc_x[,j]*sqrt(pcs_explained_var[j])}

# # All cohorts: example analysis using kmeans
# # Examine the number of clusters
# wss <- sapply(1:10,function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
# plot(1:10, wss,
#      type="b", pch = 19, frame = FALSE,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")
# set.seed(123)
# kmeans_res <- kmeans(pc_x, 5, nstart = 100)$cluster
# table(kmeans_res)
# write.table(table(kmeans_res,d$CohortName))
# table(kmeans_res,alldata_is_jap[rownames(pc_x)]) # Japanese are well clustered and removed

# hierarchical clustering
dd = dist(pc_x,method="manhattan")
h = hclust(dd,method = "single")

wss <- sapply(seq(1,5000,by=100),function(k){tot_wss_hluct(k,h,pc_x)})
plot(seq(1,5000,by=100), wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
wss <- sapply(seq(1,200,by=10),function(k){tot_wss_hluct(k,h,pc_x)})
plot(seq(1,200,by=10), wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# # Take the closest sample to each of our subjects and recalculate
# our_samples = rownames(d)[d$CohortName == "elite"]
# ukbb_samples = rownames(d)[d$CohortName == "ukbb"]
# pc_x = as.matrix(d[,paste("PC",1:10,sep="")])
# rownames(pc_x) = rownames(d)
# selected_samples = c()
# n_to_select = 1
# for (i in our_samples){
#   curr_dists = sweep(pc_x[ukbb_samples,],2,pc_x[i,])
#   curr_dists = abs(curr_dists^2)
#   curr_dists = sqrt(rowSums(curr_dists))
#   curr_dists = sort(curr_dists,decreasing = F)
#   selected_samples = union(selected_samples,names(curr_dists)[1:n_to_select])
#   print(length(selected_samples))
# }
# 
# subjects_for_analysis = c(our_samples,selected_samples)
# subjects_for_analysis = rownames(d)
# pc_ps = c()
# for(j in 1:40){
#   pc_ps[j] = compute_pc_vs_binary_variable_association_p(
#     pc = d[subjects_for_analysis,paste("PC",j,sep="")],
#     y = d[subjects_for_analysis,"CohortName"]
#   )
# }
# pc_ps = p.adjust(pc_ps)
# pc_inds = which(pc_ps < 0.01) # Before correction: almost all
# 
# inds = subjects_for_analysis

kmeans_res = run_hclust(pc_x,150,dd,h)
kmeans_res[kmeans_res!=1] = 0
table(kmeans_res)
table(kmeans_res,d$CohortName)

inds = rownames(d)
res = two_d_plot_visualize_covariate(d[inds,]$PC1,d[inds,]$PC2,kmeans_res,kmeans_res,
    main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="topleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[inds,]$PC2,d[inds,]$PC3,kmeans_res,kmeans_res,
    main = "Clustering results: PCs 2 and 3",xlab="PC3",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d[inds,]$PC1,d[inds,]$PC2,d[inds,]$CohortName,d[inds,]$CohortName,
    main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[inds,]$PC2,d[inds,]$PC3,d[inds,]$CohortName,d[inds,]$CohortName,
    main = "All cohorts",xlab="PC3",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[inds,]$PC39,d[inds,]$PC1,d[inds,]$CohortName,d[inds,]$CohortName,
    main = "All cohorts",xlab="PC3",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d[inds,]$PC39,d[inds,]$PC1,kmeans_res,kmeans_res,
                                     main = "All cohorts",xlab="PC3",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])


# res = two_d_plot_visualize_covariate(pca2[,"PC1"],pca2[,"PC2"],d[inds,]$CohortName,d[inds,]$CohortName,
#   main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
# legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
# res = two_d_plot_visualize_covariate(pca2[,"PC13"],pca2[,"PC14"],d[inds,]$CohortName,d[inds,]$CohortName,
#   main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
# legend(x="bottomleft",names(res[[1]]),fill = res[[1]])


# pcs_explained_var = read.table("merged_data_qctool_bed.eigenval")[,1]
# pcs_explained_var = pcs_explained_var/sum(pcs_explained_var)
# pcs_explained_var = format(pcs_explained_var,digits = 1)
# par(mfrow=c(3,4))
# for(j in 1:12){
#   x1 = newd[,paste("PC",j,sep="")]
#   y1 = newd$ExerciseGroup
#   currd = data.frame(x1,y1)
#   boxplot(x1~y1,data=currd,main = paste("PC",j," (",pcs_explained_var[j],")",sep=""))
#   print(paste(j,cor.test(x1,newd$ExerciseGroup)$p.value))
# }
# pc_matrix = newd[,grepl("^PC",colnames(newd))]
# vars = apply(pc_matrix,2,var)

