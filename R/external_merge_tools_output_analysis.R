
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
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/new_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/new_bed_2.frq"
our_data_mafs_by_group = list(
  "genepool" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/genepool_cohort_freq.frq",
  "cooper" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/Cooper_cohort_freq.frq",
  "elite" = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/ELITE_cohort_freq.frq"
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
names(alldata_is_jap) = d$IID
set.seed(123)
pc_x = as.matrix(d[,paste("PC",1:10,sep="")])
rownames(pc_x) = d$IID

# kmeans_res <- kmeans(pc_x, 3, nstart = 25)$cluster

run_hclust<-function(pc_x,k,dd=NULL,h=NULL){
  if(is.null(dd)){dd = dist(pc_x,method="manhattan")}
  if(is.null(h)){h = hclust(dd,method = "complete")}
  clust = cutree(h,k)
  return(clust)
}
dd = dist(pc_x,method="manhattan")
h = hclust(dd,method = "complete")
kmeans_res <- run_hclust(pc_x, 10, d=d,h=h)

table(kmeans_res)
table(kmeans_res,d[rownames(pc_x),]$CohortName)
table(kmeans_res,alldata_is_jap[rownames(pc_x)]) # Japanese are well clustered and removed

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
                 paste("--covar-name sex,Age,",paste(paste("PC",1:5,sep=""),collapse=","),sep="")
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
                 paste("--covar-name sex,Age,",paste(paste("PC",1:5,sep=""),collapse=","),sep=""),
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
                 paste("--covar-name sex,Age,",paste(paste("PC",1:5,sep=""),collapse=","),sep=""),
                 "--adjust --out",paste(out_path,"ukbb_vs_cooper_logistic",sep=''))
curr_sh_file = "ukbb_vs_cooper_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Sanity check 1 : UKBB vs. Genepool
load(paste(out_path,"clustering_data.RData",sep=""))
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
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
                 "--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--adjust --out",paste(out_path,"ukbb_vs_genepool_logistic",sep=''))
curr_sh_file = "ukbb_vs_genepool_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

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
                 "--bfile",bfile,"--logistic --flip-scan --allow-no-sex --test-missing",
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
####################################################################################################
# Analysis of the sanity check results
# mafs
get_mafs_from_file<-function(path){
  freqs = read.table(path,header=T,stringsAsFactors = F)
  rownames(freqs) = freqs[,2]
  mafs = freqs$MAF;names(mafs)=rownames(freqs)
  return(mafs)
}
combined_mafs = get_mafs_from_file(paste(bfile,".frq",sep=""))
our_mafs = get_mafs_from_file(our_data_mafs)
external_mafs = get_mafs_from_file(external_data_mafs)
names(our_mafs) = names(external_mafs)
our_group_mafs = lapply(our_data_mafs_by_group,get_mafs_from_file)
cor(our_group_mafs[[1]],our_group_mafs[[2]])

# gwas'
ukbb_vs_gp_res = read.table(
  paste(out_path,"ukbb_vs_genepool_logistic.ExerciseGroup.glm.logistic.hybrid.adjusted",sep=''),
  header=F,stringsAsFactors = F,check.names = F)
rownames(ukbb_vs_gp_res) = ukbb_vs_gp_res[,2]
colnames(ukbb_vs_gp_res) = ukbb_vs_gp_res[1,]
ukbb_vs_gp_res = ukbb_vs_gp_res[-1,]
ukbb_vs_nonukbb_res = read.table(
  paste(out_path,"ukbb_vs_nonukbb_logistic.missing.adjusted",sep=""),
  header=F,stringsAsFactors = F
)
rownames(ukbb_vs_nonukbb_res) = ukbb_vs_nonukbb_res[,2]
colnames(ukbb_vs_nonukbb_res) = ukbb_vs_nonukbb_res[1,]
ukbb_vs_nonukbb_res = ukbb_vs_nonukbb_res[-1,]
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
  snps_to_exclude_from_results=snps_to_exclude_from_results)

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
names(alldata_is_jap) = d$IID
our_pca = read_pca_res("../../analysis/maf_filter.eigenvec")

# Examine the PCA results
library(corrplot)
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
our_pca = read_pca_res("../../analysis/maf_filter.eigenvec")
all(rownames(pca1)==rownames(pca2))
corrs = cor(pca1,pca2)
corrplot(corrs)
pcainds = intersect(rownames(pca1),rownames(our_pca))
corrs = cor(pca1[pcainds,],our_pca[pcainds,])
corrplot(corrs)

run_hclust<-function(pc_x,k,dd=NULL,h=NULL){
  if(is.null(dd)){dd = dist(pc_x,method="manhattan")}
  if(is.null(h)){h = hclust(dd,method = "complete")}
  clust = cutree(h,k)
  return(clust)
}
# Cluster by PCs
# Assumption: use PC 3 and onwards as the first two are
# batch PCs that separate UKBB from non-UKBB
pc_x = as.matrix(d[,paste("PC",1:10,sep="")])
# pcs_explained_var = read.table("merged_data_qctool_bed.eigenval")[,1]
# for(j in 1:ncol(pc_x)){pc_x[,j]=pc_x[,j]*sqrt(pcs_explained_var[j])}
# apply(pc_x,2,sd)
# pc_x = pc_x[d$CohortName!="ukbb",]
dd = dist(pc_x,method="manhattan")
h = hclust(dd,method = "complete")
# # Examine the number of clusters
# wss <- sapply(1:10,function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
# wss <- sapply(2:10,function(k,d,h,pc_x){run_hclust(pc_x=pc_x, k,d=d,h=h)},pc_x,d,h)
# plot(1:10, wss,
#      type="b", pch = 19, frame = FALSE,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")

# Cluster subjects and visualize
set.seed(123)
#kmeans_res <- kmeans(pc_x, 5, nstart = 100)$cluster
kmeans_res <- run_hclust(pc_x, 20, d=d,h=h)
table(kmeans_res)
write.table(table(kmeans_res,d[rownames(pc_x),]$CohortName))
table(kmeans_res,alldata_is_jap[rownames(pc_x)]) # Japanese are well clustered and removed
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC1,
    d[rownames(pc_x),]$PC2,kmeans_res,kmeans_res,
    main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC2,
    d[rownames(pc_x),]$PC3,kmeans_res,kmeans_res,
    main = "All cohorts",xlab="PC3",ylab="PC2")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC4,
    d[rownames(pc_x),]$PC3,kmeans_res,kmeans_res,
    main = "Clustering results: PCs 3 and 4",xlab="PC3",ylab="PC4")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC4,
    d[rownames(pc_x),]$PC5,kmeans_res,kmeans_res,
    main = "All cohorts",xlab="PC4",ylab="PC5")
legend(x="top",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC6,
    d[rownames(pc_x),]$PC5,kmeans_res,kmeans_res,
    main = "Clustering results: PCs 5 and 6",xlab="PC6",ylab="PC5")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])

# Additional PCA plots
# All cohorts
table(d[rownames(pc_x),]$CohortName)
cexs = rep(1,nrow(d))
cexs[d[rownames(pc_x),]$CohortName=="ukbb"] = 0.1
cex[d[rownames(pc_x),]$CohortName == "elite"] = 1.2
inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC1[inds],
    d[rownames(pc_x),]$PC2[inds],cohorts[inds],cohorts[inds],cex=cexs,lwd=2,
    main = "All four cohors: PCs 1,2",xlab="PC1",ylab="PC2")
legend(x="bottom",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC3[inds],
    d[rownames(pc_x),]$PC2[inds],cohorts[inds],cohorts[inds],cex=cexs,lwd=2,
    main = "All four cohors: PCs 2,3",xlab="PC3",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d[rownames(pc_x),]$PC3[inds],
    d[rownames(pc_x),]$PC4[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# Check PCs vs. our cohorts in the reduced dataset
cl_tb = table(kmeans_res)
cl_fa = names(cl_tb)[cl_tb==max(cl_tb)]
selected_subjects_for_gwas = names(kmeans_res)[kmeans_res==cl_fa]
length(selected_subjects_for_gwas)
newd = d[selected_subjects_for_gwas,]
newd = d
newd = newd[newd$CohortName!="genepool",]
newd = newd[newd$CohortName!="ukbb",]
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

# direct selection of controls
our_samples = d$IID[d$CohortName=="elite" | d$CohortName == "cooper" | d$CohortName == "genepool"]
pc_x = as.matrix(our_pca[our_samples,paste("PC",1:10,sep="")])
wss <- sapply(1:20,function(k){kmeans(pc_x, k, nstart=50,iter.max = 15)$tot.withinss})
plot(1:20, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
kmeans_res = kmeans(pc_x,13)$cluster
table(kmeans_res,d[names(kmeans_res),]$CohortName)
res = two_d_plot_visualize_covariate(pc_x[,"PC1"],pc_x[,"PC2"],kmeans_res,kmeans_res,
    main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pc_x[,"PC3"],pc_x[,"PC4"],kmeans_res,kmeans_res,
    main = "Clustering results: PCs 3 and 4",xlab="PC3",ylab="PC4")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pc_x[,"PC5"],pc_x[,"PC6"],kmeans_res,kmeans_res,
    main = "Clustering results: PCs 5 and 6",xlab="PC5",ylab="PC6")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pc_x[,"PC7"],pc_x[,"PC8"],kmeans_res,kmeans_res,
    main = "Clustering results: PCs 7 and 8",xlab="PC7",ylab="PC8")
legend(x="bottomright",names(res[[1]]),fill = res[[1]])



controls = d$IID[d$CohortName=="ukbb"]

