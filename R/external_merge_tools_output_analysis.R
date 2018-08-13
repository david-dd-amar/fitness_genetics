
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
python_script = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/python_sh/recode_indels.py"
source(script_file)


bgen_file = "/oak/stanford/groups/euan/projects/ukbb/qctool/aug10_full_set"
out_path = "/oak/stanford/groups/euan/projects/ukbb/qctool/"

external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age_with_info.txt"

our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/three_group_analysis_genepool_controls.phe"
our_phe_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/three_group_analysis_genepool_controls.phe"


######################################################################################################
# Correct some issues in the sample file
sample_file = paste(bgen_file,".sample",sep="")
sample_bgen_data = read.table(sample_file,header=T,stringsAsFactors = F)
table(sample_bgen_data$sex)
table(is.na(sample_bgen_data$sex))
row_inds = which(is.na(sample_bgen_data$sex) | 
                   sample_bgen_data$sex != "1" | 
                   sample_bgen_data$sex != "2")
row_inds = setdiff(row_inds,1)
sample_bgen_data[row_inds,"sex"] = "0"
write.table(sample_bgen_data,file=sample_file,sep=" ",
            quote=F,row.names = F)

# ####################################################################################################
# The code below reruns some of the analyses using merged bgen files from qctools
# Run PCA and freq
err_path = paste(out_path,"maf_and_pca.err",sep="")
log_path = paste(out_path,"maf_and_pca.log",sep="")
curr_cmd = paste("plink2 --bgen",paste(bgen_file,".bgen",sep=""),
                 "--sample",paste(bgen_file,".sample",sep=""),
                 "--pca --freq --out",paste(out_path,"merged_bed_final_for_gwas",sep=''))
curr_sh_file = "maf_and_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000,plink_pkg = "plink/2.0a1"),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

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
genepool_inds = which(as.numeric(our_phe)==-1)
# all(our_covars[,1]==our_phe[,1])
# table(our_phe[,3])
external_covars = read.table(external_covars_path,stringsAsFactors = F)
external_covars = cbind(as.character(external_covars[,1]),external_covars)
external_samples = as.character(external_covars[,1])

# # read our fam file
# fam_info = read_plink_table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/maf_filter.fam",has_header = F)
# iid_to_fid = fam_info[,1]

# temp sol for ukbb - add batches
batch_data = as.matrix(
  read.table("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.tab.fam",
             stringsAsFactors = F, header=T))
rownames(batch_data) = as.character(batch_data[,1])
external_covars = cbind(batch_data[external_samples,c(1:2,5:6)],external_covars[,4])
external_covars = cbind(external_covars,rep("-1",nrow(external_covars)))
colnames(external_covars) = c("FID","IID","sex","Batch","Age","ExerciseGroup")

# Define the final covariance matrix (with the ExerciseGroup column)
covars = rbind(our_covars[,c(1,2,3,5,6,4)],external_covars)
pca_res = read.table(paste(out_path,"merged_bed_final_for_gwas.eigenvec",sep=''),stringsAsFactors = F)
rownames(pca_res) = pca_res[,2]
rownames(pca_res) = gsub("-","_",rownames(pca_res))
length(intersect(rownames(pca_res),covars[,"IID"]))
setdiff(rownames(pca_res),covars[,"IID"])
pca_res = pca_res[,-c(1:2)]
colnames(pca_res) = paste("PC",1:ncol(pca_res),sep='')
table(is.element(covars[,"IID"],set=rownames(pca_res))) # all should be true
covars = cbind(covars,pca_res[covars[,"IID"],])
covars = covars[-genepool_inds,]
covars[,"Batch"] = cov_phe_col_to_plink_numeric_format(covars[,"Batch"])
table(covars[,"Batch"])
table(covars[,"Age"])
for(j in 1:ncol(covars)){
  covars[,j] = gsub(" ","",as.character(covars[,j]))
}

ind = which(colnames(covars)=="ExerciseGroup")
write.table(file=paste(out_path,"ukbb_elite_cooper.phe",sep=''),
            covars[,c(1:2,ind)],sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=paste(out_path,"ukbb_elite_cooper_covars.phe",sep=''),
            covars[,-ind],sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# Run GWAS

# 1. Linear of all three groups + sex, age, and 10 PCs
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
err_path = paste(out_path,"gwas_three_groups_linear.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bgen",paste(bgen_file,".bgen",sep=""),
                 "--sample",paste(bgen_file,".sample",sep=""),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))


# 1. Linear of all three groups + sex, age, and 5 PCs
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
err_path = paste(out_path,"gwas_three_groups_linear_5pcs.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear_5pcs.log",sep="")
curr_cmd = paste("plink2",
                 "--bgen",paste(bgen_file,".bgen",sep=""),
                 "--sample",paste(bgen_file,".sample",sep=""),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear_5pcs",sep=''))
curr_sh_file = "gwas_three_groups_linear_5pcs.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################


setwd("/Users/David/Desktop/elite/gwas_results/ukbb_qctools_pca/")
d = read.table("ukbb_elite_cooper_covars.phe",header=T,stringsAsFactors = F)
d2 = read.delim("../../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
table(is_jap)
d_cohort = read.table("ukbb_elite_cooper.phe",header=T,stringsAsFactors = F)
cohorts = d_cohort$ExerciseGroup
cohorts[cohorts=="-1"] = "UKBB"
cohorts[cohorts=="1"] = "ELITE"
cohorts[cohorts=="0"] = "Cooper"
names(cohorts) = d_cohort[,2]
table(cohorts)
all(names(cohorts)==d[,2])

# Cluster by PCs
# Assumption: use PC 3 and onwards as the first two are
# batch PCs that separate UKBB from non-UKBB
pc_x = as.matrix(d[,c("PC1","PC2")])
pc_x = as.matrix(d[,c("PC3","PC2")])

# Check number of clusters in PCA plot
wss <- sapply(1:10,
              function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
wss[2:10]/wss[1:9] # by manual inspection we choose 3 clusters here (consistent with 
# the analysis of our samples without UKBB)
plot(1:10, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# PCA plots
inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d$PC1[inds],
                                     d$PC2[inds],cohorts[inds],cohorts[inds],
                                     main = "UKBB, Cooper, ELITE",xlab="PC1",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
                                     d$PC2[inds],cohorts[inds],cohorts[inds],
                                     main = "UKBB, Cooper, ELITE",xlab="PC3",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
                                     d$PC4[inds],cohorts[inds],cohorts[inds],
                                     main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

table(kmeans_res)
res = two_d_plot_visualize_covariate(d$PC12[inds],
                                     d$PC11[inds],kmeans_res,kmeans_res,
                                     main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])


