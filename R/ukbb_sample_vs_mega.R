
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Our imputed data, generated on March 2019
bfiles = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_ukbb/"
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_ukbb/new_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_ukbb/new_bed_2.frq"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_ukbb/gwas/"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"
our_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"

try({system(paste("mkdir",out_path))})
chrs=1:22

####################################################################################################
####################################################################################################
####################################################################################################
# Metadata

# Read covars, pca, and create phe file
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

# Define the merged covariance matrix (with the ExerciseGroup column), without the PCs
our_covars_wo_pcs = cbind(our_covars[,c("FID","IID","sex","batch","age","Cohort")],cohorts)
colnames(our_covars_wo_pcs) = colnames(external_covars)
covars = as.matrix(rbind(our_covars_wo_pcs,external_covars))
rownames(covars) = covars[,"IID"]
covars[covars[,7]=="2",7] = "elite"
covars[covars[,7]=="1",7] = "cooper"
covars = covars[,-6]

# columns for pairwise comparisons
gp_vs_ukbb = rep(NA,nrow(covars))
gp_vs_ukbb[covars[,"CohortName"] == "genepool"] = "1"
gp_vs_ukbb[covars[,"CohortName"] == "ukbb"] = "2"

elite_vs_ukbb = rep(NA,nrow(covars))
elite_vs_ukbb[covars[,"CohortName"] == "elite"] = "1"
elite_vs_ukbb[covars[,"CohortName"] == "ukbb"] = "2"

cooper_vs_ukbb = rep(NA,nrow(covars))
cooper_vs_ukbb[covars[,"CohortName"] == "cooper"] = "1"
cooper_vs_ukbb[covars[,"CohortName"] == "ukbb"] = "2"

covars = cbind(covars,gp_vs_ukbb,elite_vs_ukbb,cooper_vs_ukbb)
write.table(covars,file=paste(out_path,"mega_ukbb_covars_wo_pcs.phe",sep=""),
            row.names = F,col.names = T,quote=F,sep=" ")

####################################################################################################
####################################################################################################
####################################################################################################
# LD prune, PCA on all samples
ld_prune_res = paste(bfiles,"ld_prune/",sep="")
system(paste("mkdir",ld_prune_res))

# Run LD prunning on each bed file
for(chr in chrs){
  curr_name = paste("ld_",chr,sep='')
  curr_cmd = paste("plink --bfile",paste(bfiles,"chr",chr,sep=''),
                   "--threads 8",
                   "--maf 0.01",
                   "--indep-pairwise 250 10",0.5,
                   "--out",paste(ld_prune_res,"chr",chr,sep='')
  )
  run_plink_command(curr_cmd,ld_prune_res,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Extract LD-pruned beds
for(chr in chrs){
  curr_name = paste("ld_extract_",chr,sep='')
  curr_cmd = paste("plink --bfile",paste(bfiles,"chr",chr,sep=''),
                   "--threads 4",
                   "--extract", paste(ld_prune_res,"chr",chr,".prune.in",sep=''),
                   "--make-bed --out",paste(ld_prune_res,"ld_pruned_chr",chr,sep='')
  )
  run_plink_command(curr_cmd,ld_prune_res,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}

# Merge all resulting bed files into a single one
all_bfiles = paste(ld_prune_res,"ld_pruned_chr",chrs,sep='')
write.table(t(t(all_bfiles[-1])),file=paste(ld_prune_res,"all_bfiles.txt",sep=""),
            row.names=F,col.names=F,quote=F)
curr_name = paste("merge_beds",sep='')
curr_cmd = paste("plink --bfile",all_bfiles[1],
                 "--merge-list",paste(ld_prune_res,"all_bfiles.txt",sep=""),
                 "--threads 4",
                 "--make-bed --out",paste(ld_prune_res,"merged_dataset",sep='')
)
run_plink_command(curr_cmd,ld_prune_res,curr_name,
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)

# Run PCA 
curr_name = "run_pca"
curr_cmd = paste("plink --bfile",paste(ld_prune_res,"merged_dataset",sep=''),
                 "--threads 16",
                 "--pca 40",
                 "--out",paste(ld_prune_res,curr_name,sep='')
)
run_plink_command(curr_cmd,ld_prune_res,curr_name,
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)

####################################################################################################
####################################################################################################
####################################################################################################
# Load the EU manual clustering analysis, do a GWAS between GP and UKBB, adjusting for sex
# and age only
eu_samples = get(load(paste(out_path,"../manual_clustering.RData",sep="")))
fam1 = read.table(paste(out_path,"../chr1.fam",sep=""),stringsAsFactors = F)
curr_fam = fam1[is.element(fam1[,2],set=names(which(eu_samples))),1:2]
pcax = read_pca_res(paste(ld_prune_res,"run_pca.eigenvec",sep=""))
covars = read.table(paste(out_path,"mega_ukbb_covars_wo_pcs.phe",sep=""), 
                    sep = " ",header=T,stringsAsFactors = F)
rownames(covars) = covars$IID

covars = cbind(covars,pcax[rownames(covars),])
curr_cov_file = paste(out_path,"mega_ukbb_covars_with_pcs.phe",sep="")
curr_sample_file = paste(out_path,"mega_ukbb_eu_fam.fam",sep="")

# fix the covar file
covars[covars=="_"] = NA
covars$Age[is.na(as.numeric(covars$Age))] = NA

write.table(covars,file=curr_cov_file,
              row.names = F,col.names = T,quote=F,sep=" ")
write.table(curr_fam,file=curr_sample_file,
            row.names = F,col.names = F,quote=F,sep=" ")
curr_path = paste(out_path,"gp_vs_ukbb/",sep="")
system(paste("mkdir",curr_path))

# Run the GWAS: try different numbers of PCs
chrs= paste("chr",1:22,sep="")
for(num_pcs in 0:3){
  cov_string = "--covar-name sex,Age"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,Age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(bfiles,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",curr_cov_file,
                     "--pheno-name gp_vs_ukbb",
                     "--covar",curr_cov_file,
                     "--maf 0.01",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--keep", curr_sample_file,
                     "--out",paste(curr_path,"gp_vs_ukbb",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("gp_vs_ukbb",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
}

# Read the results, record the SNPs to exclude
gp_vs_ukbb_res = list()
all_files = list.files(curr_path)
all_files = all_files[grepl("adjusted",all_files)]
for(num_pcs in 0:3){
  curr_pcs_snps = c()
  for(chr in chrs){
    currf = all_files[grepl(paste(chr,'.gp',sep=""),all_files) &
                         grepl(paste("ukbb",num_pcs,sep=""),all_files)]
    print(currf)
    curr_res = read.table(paste(curr_path,currf,sep=""),stringsAsFactors = F)
    curr_snps = curr_res[curr_res[,3] <= 1e-5,2]
    print(length(curr_snps))
    curr_pcs_snps = c(curr_pcs_snps,curr_snps)
  }
  gp_vs_ukbb_res[[as.character(num_pcs)]] = curr_pcs_snps
}
sapply(gp_vs_ukbb_res,length)
save(gp_vs_ukbb_res,file=paste(out_path,"gp_vs_ukbb_res.RData",sep=""))

####################################################################################################
####################################################################################################
####################################################################################################
# Reanalyze the data using EUs only and exclude variants
# using the gp vs ukbb analysis above
# The structure of the output is as follows:
# For each option to set the number of pcs in the ukbb vs gp analysis
# we create a directory
# In this directory we will have one folder for the ld-pca analysis
# and one folder for the GWAS analyses

curr_sample_file = paste(out_path,"mega_ukbb_eu_fam.fam",sep="")
load(paste(out_path,"gp_vs_ukbb_res.RData",sep=""))

for(nn in names(gp_vs_ukbb_res)){
  
  curr_dir = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/",sep="")
  system(paste("mkdir",curr_dir))
  
  surr_snps_to_exclude = gp_vs_ukbb_res[[nn]]
  surr_snps_to_exclude_file = paste(curr_dir,"surr_snps_to_exclude.txt",sep="")
  write.table(t(t(surr_snps_to_exclude)),file=surr_snps_to_exclude_file,row.names = F,col.names = F,
              quote = F)
  
  # LD prune, PCA on all samples
  ld_prune_res = paste(curr_dir,"ld_prune/",sep="")
  system(paste("mkdir",ld_prune_res))
  
  # Run LD prunning on each bed file
  for(chr in chrs){
    curr_name = paste("ld_",chr,sep='')
    curr_cmd = paste("plink --bfile",paste(bfiles,"chr",chr,sep=''),
                     "--threads 8",
                     "--maf 0.01",
                     "--indep-pairwise 250 10",0.5,
                     "--keep", curr_sample_file,
                     "--exclude",surr_snps_to_exclude_file,
                     "--out",paste(ld_prune_res,"chr",chr,sep='')
    )
    run_plink_command(curr_cmd,ld_prune_res,curr_name,
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
  }
  wait_for_job(600)
  # Extract LD-pruned beds
  for(chr in chrs){
    curr_name = paste("ld_extract_",chr,sep='')
    curr_cmd = paste("plink --bfile",paste(bfiles,"chr",chr,sep=''),
                     "--threads 4",
                     "--keep", curr_sample_file,
                     "--exclude",surr_snps_to_exclude_file,
                     "--extract", paste(ld_prune_res,"chr",chr,".prune.in",sep=''),
                     "--make-bed --out",paste(ld_prune_res,"ld_pruned_chr",chr,sep='')
    )
    run_plink_command(curr_cmd,ld_prune_res,curr_name,
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
  }
  wait_for_job(600)
  # Merge all resulting bed files into a single one
  all_bfiles = paste(ld_prune_res,"ld_pruned_chr",chrs,sep='')
  write.table(t(t(all_bfiles[-1])),file=paste(ld_prune_res,"all_bfiles.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  curr_name = paste("merge_beds",sep='')
  curr_cmd = paste("plink --bfile",all_bfiles[1],
                   "--merge-list",paste(ld_prune_res,"all_bfiles.txt",sep=""),
                   "--threads 4",
                   "--make-bed --out",paste(ld_prune_res,"merged_dataset",sep='')
  )
  run_plink_command(curr_cmd,ld_prune_res,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
  wait_for_job(600)
  print(paste("ld prune done, number of variants for PCA:",
              nrow(read.table(paste(ld_prune_res,"merged_dataset.bim",sep='')))))
  gc()
  # Run PCA 
  curr_name = "run_pca"
  curr_cmd = paste("plink --bfile",paste(ld_prune_res,"merged_dataset",sep=''),
                   "--threads 16",
                   "--pca 40",
                   "--out",paste(ld_prune_res,curr_name,sep='')
  )
  run_plink_command(curr_cmd,ld_prune_res,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)
}


####################################################################################################
####################################################################################################
####################################################################################################

# Locally
setwd("/Users/David/Desktop/elite/feb_2019/mega_ukbb")
d2 = read.table("20k_rand_controls_sex_age_with_info.txt",stringsAsFactors=F,header=F,row.names = 1)
source("~/Desktop/repos/fitness_genetics/R/gwas_flow_helper_functions.R")
pcax = read_pca_res("run_pca.eigenvec")
d = read.table("integrated_sample_metadata_and_covariates.phe",header=T,stringsAsFactors = F)
rownames(d) = d$IID
inds = intersect(rownames(d),rownames(pcax))
d = d[inds,]
dim(d)
mega_eus = read.table("eu_pheno.phe",header=T,stringsAsFactors = F)
rownames(mega_eus) = mega_eus$IID

# Parse the cohorts
cohorts_all = rep("UKBB",nrow(pcax))
names(cohorts_all) = rownames(pcax)
cohorts_all[rownames(d)[d$Cohort=="1" & !is.na(d$Cohort)]] = "Cooper"
cohorts_all[rownames(d)[d$Cohort=="2" & !is.na(d$Cohort)]] = "ELITE"
cohorts_all[rownames(d)[d$Cohort=="genepool" & !is.na(d$Cohort)]] = "GP"
table(cohorts_all)

# PCA plots
inds = rownames(pcax)
curr_cohorts = cohorts_all
curr_cex = rep(0.2,length(inds))
curr_cex[curr_cohorts[inds]=="UKBB"] = 0.01

res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],curr_cohorts[inds],curr_cohorts[inds],
            main = "All filters",xlab="PC1",ylab="PC2",lwd=3,cex=curr_cex)
legend(x="bottomright",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],curr_cohorts[inds],curr_cohorts[inds],
            main = "All filters",xlab="PC1",ylab="PC2",lwd=3,cex=curr_cex,
            xlim=c(-0.01,0.003),ylim=c(-0.04,0.02))
legend(x="bottomright",names(res[[1]]),fill = res[[1]])

# plot a single pc
pc=1
par(mfrow=c(2,1))
hist(pcax[inds[curr_cohorts[inds]=="UKBB"],pc],breaks=200,main="PC1 UKBB samples")
hist(pcax[inds[curr_cohorts[inds]!="UKBB"],pc],breaks=200, main = "PC1 non-UKBB samples")

# Do the pcs perfectly predict the group?
y = as.factor(curr_cohorts=="ukbb")
x = pcax[,1:10]
logistic_d = data.frame(y,x)
cv_res = c()
folds = sample(rep(1:10,length(y)/10))[1:length(y)]
for(i in 1:10){
  tr = logistic_d[folds!=i,]
  te = x[folds==i,]
  tey = y[folds==i]
  m = glm(y~.,family = "binomial",data=tr)
  preds = predict(m,newdata=data.frame(logistic_d[folds==i,]),type="response")
  cv_res = rbind(cv_res,cbind(preds,tey))
}
table(cv_res[,1]>0.5,cv_res[,2])

# select EUs
manual_clustering = pcax[,1]> -0.01 & pcax[,2] > -0.02 & pcax[,2]<0.02
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],manual_clustering[inds],manual_clustering[inds],
                                     main = "All filters",xlab="PC1",ylab="PC2",lwd=3,cex=curr_cex)
table(curr_cohorts[manual_clustering])
save(manual_clustering,file="manual_clustering.RData")

# compare to our mega-only clustering
xx1 = get(load("/Users/David/Desktop/elite/november2018_analysis/mega_with_gp/manual_clustering.RData"))
xx2 = manual_clustering[names(xx1)]
table(xx1,xx2)

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
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
library(class)
pc_ps = c()
for(j in 1:40){
  # pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
  #   pc = d[subjects_for_analysis,paste("PC",j,sep="")],
  #   y = d[subjects_for_analysis,"CohortName"]
  # )
  pc_ps[j] = chisq.test(table(knn.cv(d[,paste("PC",j,sep="")],
                                     d$CohortName,k=5),d$CohortName))$p.value
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

# Elite vs. UKBB: select samples
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID

# Take the closest sample to each of our subjects and recalculate
our_samples = rownames(d)[d$CohortName == "elite"]
ukbb_samples = rownames(d)[d$CohortName == "ukbb"]
pc_x = as.matrix(d[,paste("PC",1:5,sep="")])
rownames(pc_x) = rownames(d)
selected_samples = c()
n_to_select = 10
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
pc_inds

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

####################################################################################################
####################################################################################################
####################################################################################################

# Cooper vs. UKBB: select samples
d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID

# Take the closest sample to each of our subjects and recalculate
our_samples = rownames(d)[d$CohortName == "cooper"]
ukbb_samples = rownames(d)[d$CohortName == "ukbb"]
pc_x = as.matrix(d[,paste("PC",1:5,sep="")])
rownames(pc_x) = rownames(d)
selected_samples = c()
n_to_select = 10
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
  # pc_ps[j] = compute_pc_vs_discrete_variable_association_p(
  #   pc = d[subjects_for_analysis,paste("PC",j,sep="")],
  #   y = d[subjects_for_analysis,"CohortName"]
  # )
  pc_ps[j] = chisq.test(table(knn.cv(d[subjects_for_analysis,paste("PC",j,sep="")],
                                  d$CohortName,k=5),d[subjects_for_analysis,"CohortName"]$CohortName))$p.value
}
pc_ps = p.adjust(pc_ps)
pc_inds = which(pc_ps < 0.01) # Before correction: almost all
pc_inds

# 5. Logistic: Cooper vs. UKBB, + sex, age, and 10 PCs
covars_copy = d[subjects_for_analysis,]
covar_file = paste(out_path,"ukbb_vs_cooper_matched_samples.phe",sep='')
write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_cooper_matched_samples.err",sep="")
log_path = paste(out_path,"ukbb_vs_cooper_matched_samples.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",gwas_bfile,"--logistic hide-covar firth-fallback",
                 paste("--pheno",covar_file),
                 paste("--pheno-name ExerciseGroup"),
                 paste("--covar",covar_file),
                 paste("--covar-name sex,Age,",paste(paste("PC",1:10,sep=""),collapse=","),sep=""),
                 "--adjust --out",paste(out_path,"ukbb_vs_cooper_logistic_matched_samples",sep=''))
curr_sh_file = "ukbb_vs_cooper_logistic_matched_samples.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

####################################################################################################
####################################################################################################
####################################################################################################

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
gwas_bfile = paste(out_path,"merged_data_after_pca_rl_filters",sep='')
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

####################################################################################################
####################################################################################################
####################################################################################################
# Additional analysis: QC and comparisons

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
table(m[,3]<5e-8,m[,2]<5e-8)

# Optional: discard unreliable UKBB snps
bad_ukbb_snps = read.table(
  "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/bad_ukbb_snps.txt",
  stringsAsFactors = F)
bad_ukbb_snps = bad_ukbb_snps[,1]

# Check if significant results tend to appear in blocks: find outlier problematic snps
res = read.table(paste(out_path,res_files[3],sep=''),stringsAsFactors = F)
gwas_bfile = paste(out_path,"merged_data_after_pca_rl_filters",sep='')
table(res[is.element(res[,2],set=flipscan_problematic_snps),3]<1e-8)
table(res[is.element(res[,2],set=bad_ukbb_snps),3]<1e-8)
bad_snps = union(bad_ukbb_snps,flipscan_problematic_snps)
table(res[is.element(res[,2],set=bad_snps),3]<5e-8)

create_fuma_files_for_fir(out_path,
                          paste(gwas_bfile,".bim",sep=""),
                          paste(gwas_bfile,".frq",sep=""),p = 1,maf = 0.01,
                          snps_to_exclude_from_results=bad_snps)

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

res = two_d_plot_visualize_covariate(d[inds,]$PC1,d[inds,]$PC40,d[inds,]$CohortName,d[inds,]$CohortName,
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

pcx = d[,grepl("^PC",colnames(d))]
rownames(d) = d$IID
rownames(pcx) = rownames(d)
excluded_subjects = c()
n = nrow(pcx)
while (nrow(pcx) > 0){
  curr_exclude = rep(F,nrow(pcx))
  g = d[rownames(pcx),]$CohortName
  for(j in 1:ncol(pcx)){
    curr_exclude = curr_exclude | simple_pc_bins_outlier_detection(pcx[,j],g)
    print(paste(j,sum(curr_exclude)))
  }
  if(sum(curr_exclude)==0){break}
  newx = pcx[!curr_exclude,]
  pcx = prcomp(newx,retx = T)$x
}

res = two_d_plot_visualize_covariate(pcx[,1],pcx[,2],g,g,
    main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcx[,31],pcx[,32],g,g,
      main = "Clustering results: PCs 1 and 2",xlab="PC1",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

pc_ps = c()
for(j in 1:40){
  disc_pc = cut(pcx[,j],breaks = 5)
  tb = table(g,disc_pc)
  tb = tb[,colSums(tb)>0]
  pc_ps[j] = chisq.test(tb)$p.value
}
p.adjust(pc_ps)

simple_pc_bins_outlier_detection<-function(pc,g,bins=10,pct = 0.99){
  disc_pc = cut(pc,breaks = bins)
  tb = table(g,disc_pc)
  ukbb_p = tb["ukbb",]/colSums(tb)
  bad_pc_bins = names(which(ukbb_p>pct))
  to_rem = is.element(disc_pc,set=bad_pc_bins)
  return(to_rem)
}





