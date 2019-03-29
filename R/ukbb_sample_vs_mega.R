
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Our imputed data (for mega and ukbb), generated on March 2019
bfiles = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_our_imp_ukbb/"
external_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_our_imp_ukbnew_bed_1.frq"
our_data_mafs = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_our_imp_ukbb/new_bed_2.frq"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/with_our_imp_ukbb/gwas/"
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
# GWAS: UKBB vs our cohorts

for(nn in names(gp_vs_ukbb_res)){
  curr_dir = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/",sep="")
  surr_snps_to_exclude_file = paste(curr_dir,"surr_snps_to_exclude.txt",sep="")
  # LD prune results are in this directory + PCA on the eu samples
  ld_prune_res = paste(curr_dir,"ld_prune/",sep="")
  pcax = read_pca_res(paste(ld_prune_res,"run_pca.eigenvec",sep=""))
  covars = read.table(paste(out_path,"mega_ukbb_covars_wo_pcs.phe",sep=""), 
                      sep = " ",header=T,stringsAsFactors = F)
  rownames(covars) = covars$IID
  covars = cbind(covars[rownames(pcax),],pcax)
  # Before we move on, update the current dir to a new one for the gwas results
  curr_dir = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
  system(paste("mkdir",curr_dir))
  # fix the covar file
  covars[covars=="_"] = NA
  covars$Age[is.na(as.numeric(covars$Age))] = NA
  # write to file
  curr_cov_file = paste(curr_dir,"mega_ukbb_covars_with_pcs.phe",sep="")
  write.table(covars,file=curr_cov_file,
              row.names = F,col.names = T,quote=F,sep=" ")
  # Run the GWAS: try different numbers of PCs
  chrs= paste("chr",1:22,sep="")
  for(num_pcs in 4:5){
    cov_string = paste("--covar-name sex,Age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")

    # elite
    for(chr in chrs){
      curr_cmd = paste("plink2 --bfile",paste(bfiles,chr,sep=''),
                       "--logistic firth-fallback hide-covar",
                       "--pheno",curr_cov_file,
                       "--pheno-name elite_vs_ukbb",
                       "--covar",curr_cov_file,
                       "--maf 0.01",
                       cov_string,
                       "--allow-no-sex --adjust",
                       "--keep", curr_sample_file,
                       "--exclude",surr_snps_to_exclude_file,
                       "--threads",4,
                       "--out",paste(curr_dir,"elite_vs_ukbb_pcs",num_pcs,"_",chr,sep="")
      )
      run_plink_command(curr_cmd,curr_dir,paste("elite_vs_ukbb",num_pcs,"_",chr,sep=""),
                        get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
    
      
    }
  }
}


####################################################################################################
####################################################################################################
####################################################################################################

# Locally
setwd("/Users/David/Desktop/elite/feb_2019/mega_ukbb")
d2 = read.table("20k_rand_controls_sex_age_with_info.txt",stringsAsFactors=F,header=F,row.names = 1)
source("~/Desktop/repos/fitness_genetics/R/gwas_flow_helper_functions.R")

# pca on all
pcax = read_pca_res("run_pca.eigenvec")
# pca on eu only, after gp_vs_ukbb_filter
pcax = read_pca_res("../mega_ukbb_eu/run_pca.eigenvec")

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





