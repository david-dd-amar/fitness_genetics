
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
python_script = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/python_sh/recode_indels.py"
source(script_file)

# assumption: merged bed has frq and pca results
bfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/merged_data_qctool_bed"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"
# this should have our original pca results: important for us
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/three_group_analysis_genepool_controls.phe"
our_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/gwas/"
try({system(paste("mkdir",out_path))})

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
our_phe = as.character(our_covars[,4])
table(our_covars[,4])
cohorts = rep("genepool",nrow(our_covars))
cohorts[our_covars[,4]=="0"] = "cooper"
cohorts[our_covars[,4]=="1"] = "elite"
our_covars[,4] = as.numeric(our_covars[,4]) + 2

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
our_covars_wo_pcs = cbind(our_covars[,c(1,2,3,5,6,4)],cohorts)
colnames(our_covars_wo_pcs) = colnames(external_covars)
covars = rbind(our_covars_wo_pcs,external_covars)
rownames(covars) = covars[,"IID"]

# Define the PCAs: the combined dataset and ours
read_pca_res<-function(path){
  pca1 = read.table(path)
  r = pca1[,2]
  rownames(pca1)=r
  pca1 = pca1[,-c(1:2)]
  pca1 = as.matrix(pca1)
  if(!mode(pca1)=="numeric"){
    print("ERROR: pca matrix is not numeric");return(NULL)
  }
  colnames(pca1) = paste("PC",1:ncol(pca1),sep="")
  return(pca1)
}
combined_pcs = read_pca_res(paste(bfile,".eigenvec",sep=""))
our_dataset_pcs = our_covars[,grepl("^PC",colnames(our_covars))]
rownames(our_dataset_pcs) = our_covars[,"IID"]
pcs_explained_var = read.table(paste(bfile,".eigenval",sep=""))

subjects_for_analysis = intersect(covars[,"IID"],rownames(combined_pcs))
covars = cbind(covars[subjects_for_analysis,],combined_pcs[subjects_for_analysis,])
covars[,"Batch"] = cov_phe_col_to_plink_numeric_format(covars[,"Batch"])
for(j in 1:ncol(covars)){
  covars[,j] = gsub(" ","",as.character(covars[,j]))
}

write.table(file=paste(out_path,"all_cohorts.phe",sep=''),
            covars,sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
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
pc_x = as.matrix(d[,c("PC1","PC2","PC3")])
rownames(pc_x) = d$IID
set.seed(123)
kmeans_res <- kmeans(pc_x, 3, nstart = 25)$cluster
table(kmeans_res)
table(kmeans_res,d$CohortName)
table(kmeans_res,alldata_is_jap) # Japanese are well clustered and removed
# get the largest cluster and take its subjects
cl_tb = table(kmeans_res)
cl_fa = names(cl_tb)[cl_tb==max(cl_tb)]
selected_subjects_for_gwas = names(kmeans_res)[kmeans_res==cl_fa]
d = d[selected_subjects_for_gwas,]
d = d[d$CohortName!="genepool",]
dim(d)
save(selected_subjects_for_gwas,kmeans_res,pc_x,file=paste(out_path,"clustering_data.RData",sep=""))
write.table(file=paste(out_path,"kmeans_cleaned.phe",sep=''),
            d,sep=" ",row.names = F,col.names = T,quote=F)
####################################################################################################
####################################################################################################
####################################################################################################
# Run GWAS
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
# 1. Linear of all three groups + sex, age, and 5 PCs
covar_file = paste(out_path,"kmeans_cleaned.phe",sep='')
err_path = paste(out_path,"gwas_three_groups_linear.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",bfile,"--linear hide-covar",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 "--covar-name sex,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 2. Logistic: Elite vs. UKBB, + sex, age, and 5 PCs
covars_copy = d[d$CohortName!="cooper",]
covars_copy$ExerciseGroup[covars_copy$ExerciseGroup=="3"] = 2
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_elite.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_elite.err",sep="")
log_path = paste(out_path,"ukbb_vs_elite.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 "--covar-name sex,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust --out",paste(out_path,"ukbb_vs_elite_logistic",sep=''))
curr_sh_file = "ukbb_vs_elite_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 3. Logistic: Cooper vs. UKBB, + sex, age, and 5 PCs
covars_copy = d[d$CohortName!="elite",]
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_cooper.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_cooper.err",sep="")
log_path = paste(out_path,"ukbb_vs_cooper.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 "--covar-name sex,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust --out",paste(out_path,"ukbb_vs_cooper_logistic",sep=''))
curr_sh_file = "ukbb_vs_cooper_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Sanity check 1 : UKBB vs. Genepool
load(paste(out_path,"clustering_data.RData",sep=""))
d = read.table(paste(out_path,"all_cohorts.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
d = d[selected_subjects_for_gwas,]
covars_copy = d[d$CohortName!="elite" & d$CohortName!="cooper",]
covars_copy$ExerciseGroup[covars_copy$CohortName=="ukbb"]="1"
covars_copy$ExerciseGroup[covars_copy$CohortName=="genepool"]="2"
table(covars_copy$ExerciseGroup)
covar_file = paste(out_path,"ukbb_vs_genepool.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_genepool.err",sep="")
log_path = paste(out_path,"ukbb_vs_genepool.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 "--covar-name sex,PC1,PC2,PC3,PC4,PC5",
                 "--adjust --out",paste(out_path,"ukbb_vs_genepool_logistic",sep=''))
curr_sh_file = "ukbb_vs_genepool_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Sanity check 2: UKBB vs. not UKBB: flip scan
d = read.table(paste(out_path,"all_cohorts.phe",sep=''),header=T,stringsAsFactors = F)
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
                 "--bfile",bfile,"--logistic hide-covar firth-fallback --flip-scan",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 "--covar-name sex,PC1,PC2,PC3,PC4,PC5",
                 "--adjust --out",paste(out_path,"ukbb_vs_nonukbb_logistic",sep=''))
curr_sh_file = "ukbb_vs_nonukbb_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size = 10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Print all GWAS results in FUMA's format

from_our_sol_to_fuma_res<-function(assoc_file,bim_file,freq_file=NULL,maf = 0.001,p=1){
  res = read.delim(assoc_file,stringsAsFactors = F)
  mafs = read.table(freq_file,stringsAsFactors = F,header=F)
  bim = read.delim(bim_file,stringsAsFactors = F,header=F)
  rownames(bim) = bim[,2]
  rownames(mafs) = mafs[,2]
  rownames(res) = res$ID
  selected_snps = intersect(
    rownames(res)[res$UNADJ <= p],
    rownames(mafs)[mafs[,5] >= maf]
  )
  m = cbind(as.character(bim[selected_snps,1]),as.character(bim[selected_snps,4]),res[selected_snps,]$UNADJ)
  colnames(m) = c("chromosome","position","P-value")
  return(m)
}

create_fuma_files_for_fir(out_path,
  paste(out_path,"merged_bed_final_for_gwas.bim",sep=""),
  paste(out_path,"merged_bed_final_for_gwas.afreq",sep=""))

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

setwd("/Users/David/Desktop/elite/gwas_results/ukbb_qctools_pca/")
d = read.table("all_cohorts.phe",header=T,stringsAsFactors = F)
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

read_pca_res<-function(path){
  pca1 = read.table(path)
  r = pca1[,2]
  rownames(pca1)=r
  pca1 = pca1[,-c(1:2)]
  pca1 = as.matrix(pca1)
  if(!mode(pca1)=="numeric"){
    print("ERROR: pca matrix is not numeric");return(NULL)
  }
  colnames(pca1) = paste("PC",1:ncol(pca1),sep="")
  return(pca1)
}
pca1 = read_pca_res("merged_data_plink.eigenvec")
pca2 = read_pca_res("merged_data_qctool_bed.eigenvec")
pca2 = pca2[rownames(pca1),]
all(rownames(pca1)==rownames(pca2))
#pca1 = pca1[d[,2],]
corrs = cor(pca1,pca2)
library(corrplot)
corrplot(corrs)
diag(corrs)

# Cluster by PCs
# Assumption: use PC 3 and onwards as the first two are
# batch PCs that separate UKBB from non-UKBB
pc_x = as.matrix(d[,c("PC1","PC2","PC3")])
# Cluster subjects and visualize
set.seed(123)
kmeans_res <- kmeans(pc_x, 3, nstart = 25)$cluster
table(kmeans_res)
table(kmeans_res,d$CohortName)
table(kmeans_res,alldata_is_jap) # Japanese are well clustered and removed
res = two_d_plot_visualize_covariate(d$PC2[inds],
                                     d$PC1[inds],kmeans_res,kmeans_res,
                                     main = "All cohorts",xlab="PC2",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC2[inds],
                                     d$PC3[inds],kmeans_res,kmeans_res,
                                     main = "All cohorts",xlab="PC3",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

# Additional PCA plots
inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC1",ylab="PC2")
legend(x="bottom",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC2[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC3",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC4[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# Check PCs vs. our cohorts in the reduced dataset
cl_tb = table(kmeans_res)
cl_fa = names(cl_tb)[cl_tb==max(cl_tb)]
selected_subjects_for_gwas = names(kmeans_res)[kmeans_res==cl_fa]
newd = d[selected_subjects_for_gwas,]
newd = newd[newd$CohortName!="genepool",]
dim(newd)
sapply(d,mode)
# check associations between PCs and exercise
pcs_explained_var = read.table("merged_data_qctool_bed.eigenval")[,1]
pcs_explained_var = pcs_explained_var/sum(pcs_explained_var)
pcs_explained_var = format(pcs_explained_var,digits = 1)
par(mfrow=c(3,4))
for(j in 1:12){
  x1 = newd[,paste("PC",j,sep="")]
  y1 = newd$ExerciseGroup
  currd = data.frame(x1,y1)
  boxplot(x1~y1,data=currd,main = paste("PC",j," (",pcs_explained_var[j],")",sep=""))
  print(paste(j,cor.test(x1,newd$ExerciseGroup)$p.value))
}
pc_matrix = newd[,grepl("^PC",colnames(newd))]
vars = apply(pc_matrix,2,var)


