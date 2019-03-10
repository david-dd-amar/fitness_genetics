# This script implements what we call "QC1" in our flow (see the presentation).
# The goal here is to perform the safest analysis possible when comparing
# our MEGA data to external databases (1000G and UKBB).
#
# The current analysis is semi-manual, comments within the code are very important.
# In the future, if needed, we shall transform this into an automatic analysis.
#
# Briefly, we do a quick QC analysis. We match datasets by looking at
# SNP ids only, we exclude SNPs that may require flipping and also palindromic SNPs.
# We then do population analysis and GWAS.
#
# ASSUMPTIONS: dataset1 is 1000g, dataset2 is mega-related with or without running
# check_bim. 

# Load auxiliary functions
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Define the datasets to be analyzed: vs 1000g
dataset1 = "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/"
dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/by_chr/"
chrs = paste("chr",c(1:22),sep="")
# define the output path
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/qc/1000g_vs_direct_mega/"
system(paste("mkdir",out_path))
external_data_cov_file = "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/integrated_call_samples_v3.20130502.ALL.panel"

# Define the datasets to be analyzed: vs ukbb
dataset1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_direct_20k_rand_controls_sex_age/"
dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/by_chr/"
chrs = paste("chr",c(1:22),sep="")
# define the output path
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/qc/ukbb_vs_direct_mega/"
system(paste("mkdir",out_path))
external_data_cov_file = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"

mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"

is_snp_paly<-function(x){
	return(x=="AT" || x=="TA" || x=="CG" || x=="GC")
}

# Analyze the datasets and prepare for a merge
for(chr in chrs){
	bim1 = read.table(paste(dataset1,chr,".bim",sep=""),stringsAsFactors=F)
	bim2 = read.table(paste(dataset2,chr,".bim",sep=""),stringsAsFactors=F)
	raw_id_intersect = intersect(bim1[,2],bim2[,2])
	bim1_dups = names(which(table(bim1[,2])>1))
	bim2_dups = names(which(table(bim2[,2])>1))
	raw_id_intersect = setdiff(raw_id_intersect,bim1_dups)
	raw_id_intersect = setdiff(raw_id_intersect,bim2_dups)
	bim1 = bim1[is.element(bim1[,2],set=raw_id_intersect),]
	bim2 = bim2[is.element(bim2[,2],set=raw_id_intersect),]
	rownames(bim1) = bim1[,2]
	rownames(bim2) = bim2[,2]
	bim1 = bim1[raw_id_intersect,]
	bim2 = bim2[raw_id_intersect,]
	all(bim1[,2]==bim2[,2]) # sanity check
	
	# Exclude palindromic SNPs
	alleles2 = paste(bim2[,5],bim2[,6],sep="")
	pali_snps = raw_id_intersect[sapply(alleles2,is_snp_paly)]
	raw_id_intersect = setdiff(raw_id_intersect,pali_snps)
	bim1 = bim1[raw_id_intersect,]
	bim2 = bim2[raw_id_intersect,]
	
	alleles1 = paste(bim1[,5],bim1[,6],sep="")
	rev_alleles1 = paste(bim1[,6],bim1[,5],sep="")
	alleles2 = paste(bim2[,5],bim2[,6],sep="")
	table(rev_alleles1 == alleles2) # almost all should be TRUE
	
	# NOTE: if we want to add the ability to handle snp flip for the merge
	# then we need another file with a list of snps to flip in one of the files
	# this will also include matching the alleles to make sure we get the correct
	# ones for the flip
	snps_to_take = raw_id_intersect[rev_alleles1 == alleles2 | alleles1 == alleles2]
	curr_snps_file = paste(out_path,chr,"_curr_snps.txt",sep="")
	write.table(t(t(snps_to_take)),file=curr_snps_file,row.names=F,col.names=F,quote=F)
	print(paste("chromosome:",chr,"number of snps:",length(snps_to_take)))
	# create a bed from dataset1
	curr_name = paste("dataset1_",chr,"_make_bed",sep='')
  	curr_cmd = paste("plink --bfile",paste(dataset1,chr,sep=""),
                   "--extract",curr_snps_file,
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"dataset1_",chr,sep='')
  	)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
    
    # create a bed from dataset2
    curr_name = paste("dataset2_",chr,"_make_bed",sep='')
  	curr_cmd = paste("plink --bfile",paste(dataset2,chr,sep=""),
                   "--extract",curr_snps_file,
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"dataset2_",chr,sep='')
  	)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(240)

# Compute freqs
for(chr in chrs){
	f1 = paste(out_path,"dataset1_",chr,sep='')
	f2 = paste(out_path,"dataset2_",chr,sep='')
	
	curr_name = paste("dataset1_",chr,"_frq",sep='')
  	curr_cmd = paste("plink --bfile",f1,
                   "--freq","--threads 4",
                   "--out",f1)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
    
	curr_name = paste("dataset2_",chr,"_frq",sep='')
  	curr_cmd = paste("plink --bfile",f2,
                   "--freq","--threads 4",
                   "--out",f2)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)	
}

# Do we have SNPs that are essentially zero maf in one dataset and 
# high in the other?
maf_thr1 = 0.001;maf_thr2 = 0.05
problematic_snps = c()
for(chr in chrs){
  f1 = paste(out_path,"dataset1_",chr,".frq",sep='')
  f2 = paste(out_path,"dataset2_",chr,".frq",sep='')
  freqs1 = read.table(f1,header=T)
  freqs2 = read.table(f2,header=T)
  rownames(freqs1) = freqs1$SNP
  rownames(freqs2) = freqs2$SNP
  inds = intersect(freqs1$SNP,freqs2$SNP)
  x1 = freqs1[inds,"MAF"]
  x2 = freqs2[inds,"MAF"]
  problematic_snps = c(problematic_snps,
                       inds[(x1 >=maf_thr2 & x2 < maf_thr1) | (x1 < maf_thr1 & x2 >= maf_thr2)]
  )
}
length(problematic_snps)

# Merge each chromosome data
for(chr in chrs){
	f1 = paste(out_path,"dataset1_",chr,sep='')
	f2 = paste(out_path,"dataset2_",chr,sep='')
	curr_name = paste("merge_",chr,sep='')
	curr_cmd = paste("plink --bfile",f1,
                   "--bmerge",f2,"--threads 8",
                   "--make-bed --out",paste(out_path,"merged_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Run LD prunning on each new file
for(chr in chrs){
	curr_name = paste("ld_",chr,sep='')
	curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--threads 8",
                   "--maf 0.01",
                   "--indep-pairwise 250 10",0.5,
                   "--out",paste(out_path,"merged_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Compare pruned results to the problematic SNPs above
all_selected_snps_for_pca = c()
for(chr in chrs){
  curr_snps = read.table(paste(out_path,"merged_",chr,".prune.in",sep=''),stringsAsFactors = F)[,1]
  all_selected_snps_for_pca = c(all_selected_snps_for_pca,curr_snps)
}
print(length(all_selected_snps_for_pca))
print(intersect(all_selected_snps_for_pca,problematic_snps))

# Extract LD-pruned beds
for(chr in chrs){
	curr_name = paste("ld_extract_",chr,sep='')
	curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--threads 8",
                   "--extract", paste(out_path,"merged_",chr,".prune.in",sep=''),
                   "--make-bed --out",paste(out_path,"merged_ld_pruned_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Merge all resulting bed files into a single one
all_bfiles = paste(out_path,"merged_ld_pruned_",chrs,sep='')
write.table(t(t(all_bfiles[-1])),file=paste(out_path,"all_bfiles.txt",sep=""),
	row.names=F,col.names=F,quote=F)
curr_name = paste("merge_beds",sep='')
curr_cmd = paste("plink --bfile",all_bfiles[1],
                   "--merge-list",paste(out_path,"all_bfiles.txt",sep=""),
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"merged_dataset",sep='')
)
run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)

# Run PCA 
curr_name = "run_pca"
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
        	    "--threads 16",
                "--pca 40",
                "--out",paste(out_path,curr_name,sep='')
)
run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)

# ADMIXTURE
adm_path = paste(out_path,"admixture_output","/",sep="")
system(paste("mkdir",adm_path))
setwd(adm_path)
for (k in 1:10){
  # 1. Run the algorithm on a number of Ks
  err_path = paste("run_admix_",k,".err",sep="")
  log_path = paste("run_admix_",k,".log",sep="")
  curr_cmd = paste(
    "/home/users/davidama/apps/admixture/admixture_linux-1.3.0/admixture",
    paste(out_path,"merged_dataset.bed",sep=''), k,"--cv",
    "-s 123 -j8 -C 0.1"
  )
  curr_sh_file = paste("run_",k,".sh",sep="")
  print_sh_file(curr_sh_file,
                get_sh_prefix_one_node_specify_cpu_and_mem(
                  err_path,log_path,Ncpu=8,mem_size=32000,time="36:00:00"),curr_cmd)
  system(paste("sbatch",curr_sh_file))  
}

# #####################################################################################
# #####################################################################################
# #####################################################################################
# # Run GWAS using selected admixture results 
# selected_admx_result = paste(out_path,"admixture_output/merged_dataset.8.Q",sep="")
# admx = read.table(selected_admx_result)
# colnames(admx) = paste("admx",1:ncol(admx),sep="")
# fam = read.table(paste(out_path,"merged_dataset.fam",sep=""),stringsAsFactors = F,row.names = 2)
# rownames(admx) = rownames(fam)
# 
# # just for comparison, look at the pca
# pcax = read_pca_res(paste(out_path,"run_pca.eigenvec",sep=""))
# nrow(admx) == nrow(pcax)
# all(rownames(admx)==rownames(pcax))
# cor(pcax[,1:5],admx[,1:5])
# 
# admx = cbind(fam[,1:2],admx)
# d2 = read.table(
#   "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/integrated_call_samples_v3.20130502.ALL.panel"
#   ,stringsAsFactors=F,header=T,row.names = 1)
# mega_covars = read.table(mega_covars_path,header=T)
# rownames(mega_covars) = mega_covars[,2]
# 
# # put the covars in a single file
# sex = rep(NA,nrow(admx));names(sex) = rownames(admx)
# inds1 = intersect(names(sex),rownames(mega_covars))
# sex[inds1] = mega_covars[inds1,"sex"]
# inds2 = intersect(rownames(d2),names(sex))
# sex[inds2] = d2[inds2,3]
# sex[sex=="female"]="2";sex[sex=="male"]="1"
# table(sex)
# 
# curr_cohorts = c(rownames(mega_covars),rownames(d2))
# names(curr_cohorts) = curr_cohorts
# curr_cohorts[rownames(mega_covars)] = mega_covars$Cohort
# curr_cohorts[rownames(d2)] = d2[,1]
# # a binary column for cooper vs 1000g
# cooper_col = rep(NA,nrow(admx))
# names(cooper_col) = rownames(admx)
# cooper_col[curr_cohorts=="1"]  = 1
# cooper_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
# # a binary column for elite vs 1000g
# elite_col = rep(NA,nrow(admx))
# names(elite_col) = rownames(admx)
# elite_col[curr_cohorts=="2"]  = 1
# elite_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
# # a binary column for gp vs 1000g
# gp_col = rep(NA,nrow(admx))
# names(gp_col) = rownames(gp_col)
# gp_col[curr_cohorts=="3"]  = 1
# gp_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
# 
# covs = cbind(admx,sex,cooper_col,elite_col,gp_col)
# colnames(covs)[1:2] = c("FID","IID")
# covs[,2] = rownames(covs)
# covs_file = paste(curr_path,"admx_subjects_covs.phe",sep="")
# curr_path = paste(out_path,"admixture_gwas/",sep="")
# system(paste("mkdir",curr_path))
# write.table(covs,file=covs_file,
#             row.names = F,col.names = T,quote=F,sep=" ")
# 
# cov_string = "--covar-name sex,admx1,admx2,admx3,admx4,admx5,admx6,admx7"
# for(chr in chrs){
#   curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
#                    "--logistic hide-covar",
#                    "--pheno",covs_file,
#                    "--pheno-name cooper_col",
#                    "--covar",covs_file,
#                    "--maf 0.01",
#                    cov_string,
#                    "--allow-no-sex --adjust",
#                    "--threads",4,
#                    "--out",paste(curr_path,"cooper_gwas_res",chr,sep="")
#   )
#   run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res",chr,sep=""),
#                     get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
# }
# for(chr in chrs){
#   curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
#                    "--logistic hide-covar",
#                    "--pheno",covs_file,
#                    "--pheno-name elite_col",
#                    "--covar",covs_file,
#                    "--maf 0.01",
#                    cov_string,
#                    "--allow-no-sex --adjust",
#                    "--threads",4,
#                    "--out",paste(curr_path,"elite_gwas_res",chr,sep="")
#   )
#   run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res",chr,sep=""),
#                     get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
# }
# for(chr in chrs){
#   curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
#                    "--logistic hide-covar",
#                    "--pheno",covs_file,
#                    "--pheno-name gp_col",
#                    "--covar",covs_file,
#                    "--maf 0.01",
#                    cov_string,
#                    "--allow-no-sex --adjust",
#                    "--threads",4,
#                    "--out",paste(curr_path,"gp_gwas_res",chr,sep="")
#   )
#   run_plink_command(curr_cmd,curr_path,paste("gp_gwas_res",chr,sep=""),
#                     get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
# }
# wait_for_job(240)
# 
# concatenate_res_files<-function(res_files,res_file){
#   for(j in 1:length(res_files)){
#     if(j==1){
#       system(paste("less",res_files[j],">",res_file))
#     }
#     if(j>1){
#       system(paste("less",res_files[j],"| grep -v SNP >>",res_file))
#     }
#   }
# }
# 
# get_lambda_from_log_file<-function(f){
#   l = readLines(f)
#   s = l[grepl(" = ",l)][1]
#   if(length(s)==0 || is.na(s) || nchar(s)==0){return(NA)}
#   s = strsplit(s,split=" = ")[[1]][2]
#   s = substr(s,start = 1,stop = 5)
#   return(as.numeric(s))
# }
# 
# res_files = paste(curr_path,"elite_gwas_res",chrs,".assoc.logistic",sep="")
# res_file = paste(curr_path,"elite_gwas_res_all.assoc",sep="")
# concatenate_res_files(res_files,res_file)
# res_files = paste(curr_path,"cooper_gwas_res",chrs,".assoc.logistic",sep="")
# res_file = paste(curr_path,"cooper_gwas_res_all.assoc",sep="")
# concatenate_res_files(res_files,res_file)
# res_files = paste(curr_path,"gp_gwas_res",chrs,".assoc.logistic",sep="")
# res_file = paste(curr_path,"gp_gwas_res_all.assoc",sep="")
# concatenate_res_files(res_files,res_file)
# 
# all_log_files = list.files(curr_path)
# all_log_files = all_log_files[grepl("log$",all_log_files)]
# setwd(curr_path)
# all_lambdas = sapply(all_log_files,get_lambda_from_log_file)
# for(cc in c("elite","cooper","gp")){
#   curr_files = all_log_files[grepl(cc,all_log_files)]
#   print(paste(cc,mean(all_lambdas[curr_files],na.rm=T)))
# }
# 
# file2pvals = c()
# res_file = paste(curr_path,"elite_gwas_res_all.assoc",sep="")
# res_file2 = paste(curr_path,"fuma_elite_gwas_res_all.assoc",sep="")
# res = read.table(res_file,header=T,stringsAsFactors = F)
# res = res[,c("SNP","P")]
# colnames(res) = c("rsID","P-value")
# write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
# p = res[,2];names(p) = res[,1]
# file2pvals[[res_file]] = p
# res_file = paste(curr_path,"cooper_gwas_res_all.assoc",sep="")
# res_file2 = paste(curr_path,"fuma_cooper_gwas_res_all.assoc",sep="")
# res = read.table(res_file,header=T,stringsAsFactors = F)
# res = res[,c("SNP","P")]
# colnames(res) = c("rsID","P-value")
# write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
# p = res[,2];names(p) = res[,1]
# file2pvals[[res_file]] = p
# res_file = paste(curr_path,"gp_gwas_res_all.assoc",sep="")
# res_file2 = paste(curr_path,"fuma_gp_gwas_res_all.assoc",sep="")
# res = read.table(res_file,header=T,stringsAsFactors = F)
# res = res[,c("SNP","P")]
# colnames(res) = c("rsID","P-value")
# write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
# p = res[,2];names(p) = res[,1]
# file2pvals[[res_file]] = p
# 
# table_after_name_matching <-function(x1,x2,thr){
#   inds = intersect(names(x1),names(x2))
#   table(x1[inds] < thr, x2[inds]<thr)
# }
# table_after_name_matching(file2pvals[[3]],file2pvals[[2]],1e-6)

#####################################################################################
#####################################################################################
#####################################################################################
curr_path = paste(out_path,"gwas_eu/",sep="")

pcax = read_pca_res(paste(out_path,"run_pca.eigenvec",sep=""))
# For 1000G metadata
d2 = read.table(external_data_cov_file,stringsAsFactors=F,header=T,row.names = 1)
# For UKBB metadata
d2 = read.table(external_data_cov_file,stringsAsFactors=F,header=F,row.names = 1)

mega_covars = read.table(mega_covars_path,header=T)
rownames(mega_covars) = mega_covars[,2]

# Define pheno: 1000G
cohorts2 = c(rownames(mega_covars),d2[,1])
names(cohorts2) = cohorts2
cohorts2[rownames(mega_covars)] = mega_covars$Cohort
cohorts2[rownames(d2)] = d2[,1]

# Define EUR-descendants manually - 1000G QC1
yy = pcax[,1] < -0.004 & pcax[,2]< -0.003
table(yy)
table(yy,cohorts2[names(yy)])

# Define pheno: UKBB
cohorts2 = c(rownames(mega_covars),rownames(d2))
names(cohorts2) = cohorts2
cohorts2[rownames(mega_covars)] = as.character(mega_covars$Cohort)
cohorts2[rownames(d2)] = "ukbb"
cohorts2[cohorts2=="1"] = "cooper"
cohorts2[cohorts2=="2"] = "elite"
table(cohorts2)

# Define EUR-descendants manually - UKBB QC1
yy = pcax[,1] < 0.005 & pcax[,2]< 0.03 & pcax[,2] > -0.03
table(yy,cohorts2[names(yy)])

# Define the selected subjects for PCA rerun
our_subjects = names(yy)[yy]
curr_fam = read.table(paste(out_path,"merged_dataset.fam",sep=''))
curr_fam = curr_fam[is.element(curr_fam[,2],set=our_subjects),]
system(paste("mkdir",curr_path))
write.table(curr_fam,file=paste(curr_path,"our_subjects.txt",sep=""),
            row.names = F,col.names = F,quote=F,sep=" ")

# Rerun the PCA, add curr maf filter
# Run PCA 
curr_name = "rerun_pca"
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
                 "--threads 16",
                 "--keep",paste(curr_path,"our_subjects.txt",sep=""),
                 "--maf 0.05",
                 "--pca 40",
                 "--out",paste(curr_path,curr_name,sep='')
)
run_plink_command(curr_cmd,curr_path,curr_name,
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)
wait_for_job()

# Put the covars in a single file: 1000G
new_pcax  = read_pca_res(paste(curr_path,"rerun_pca.eigenvec",sep=""))
rownames(curr_fam) = curr_fam[,2]
sex = rep(NA,nrow(new_pcax))
names(sex) = rownames(new_pcax)
inds1 = intersect(names(sex),rownames(mega_covars))
sex[inds1] = mega_covars[inds1,"sex"]
inds2 = intersect(rownames(d2),names(sex))
sex[inds2] = d2[inds2,3]
table(sex)
sex[sex=="female"]="2"
sex[sex=="male"]="1"
table(sex)
# Define the subcohorts and the pheno files
curr_cohorts = cohorts2[rownames(new_pcax)]
cooper_col = rep(NA,nrow(new_pcax))
names(cooper_col) = rownames(new_pcax)
cooper_col[curr_cohorts=="1"]  = 1
cooper_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
elite_col = rep(NA,nrow(new_pcax))
names(elite_col) = rownames(new_pcax)
elite_col[curr_cohorts=="2"]  = 1
elite_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
gp_col = rep(NA,nrow(new_pcax))
names(gp_col) = rownames(new_pcax)
gp_col[curr_cohorts=="3"]  = 1
gp_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2

# Put the covars in a single file: UKBB
# Jan 2019: Notice that "2" in our UKBB file is male - this does not fit plink's code
# so we reverse it in the script below
new_pcax  = read_pca_res(paste(curr_path,"rerun_pca.eigenvec",sep=""))
rownames(curr_fam) = curr_fam[,2]
sex = rep(NA,nrow(new_pcax))
names(sex) = rownames(new_pcax)
inds1 = intersect(names(sex),rownames(mega_covars))
sex[inds1] = mega_covars[inds1,"sex"]
inds2 = intersect(rownames(d2),names(sex))
sex[inds2] = d2[inds2,1]
table(sex)
sex[sex=="female"]="1"
sex[sex=="male"]="2"
table(sex)
# Define the subcohorts and the pheno files
curr_cohorts = cohorts2[rownames(new_pcax)]
cooper_col = rep(NA,nrow(new_pcax))
names(cooper_col) = rownames(new_pcax)
cooper_col[curr_cohorts == "cooper"]  = 1
cooper_col[curr_cohorts == "ukbb"]  = 2
elite_col = rep(NA,nrow(new_pcax))
names(elite_col) = rownames(new_pcax)
elite_col[curr_cohorts=="elite"]  = 1
elite_col[curr_cohorts == "ukbb"]  = 2
gp_col = rep(NA,nrow(new_pcax))
names(gp_col) = rownames(new_pcax)
gp_col[curr_cohorts=="genepool"]  = 1
gp_col[curr_cohorts == "ukbb"]  = 2

covs = cbind(curr_fam[rownames(new_pcax),],sex,new_pcax,cooper_col,elite_col,gp_col)
colnames(covs)[1:2] = c("FID","IID")
covs_file = "our_subjects_covs.phe"
write.table(covs,file=paste(curr_path,covs_file,sep=""),
            row.names = F,col.names = T,quote=F,sep=" ")

curr_path = paste(out_path,"gwas_eu/",sep="")
covs_file = paste(curr_path,"our_subjects_covs.phe",sep="")

# Run the GWAS: try different numbers of PCs
for(num_pcs in 0:5){
  cov_string = "--covar-name sex"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  for(chr in chrs){
    curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                     "--logistic hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name cooper_col",
                     "--covar",covs_file,
                     "--maf 0.01",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
  }
  for(chr in chrs){
    curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                     "--logistic hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name elite_col",
                     "--covar",covs_file,
                     "--maf 0.01",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
  }
  for(chr in chrs){
    curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                     "--logistic hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name gp_col",
                     "--covar",covs_file,
                     "--maf 0.01",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"gp_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("gp_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
  }
}

# Merge the results into single files
concatenate_res_files<-function(res_files,res_file){
  for(j in 1:length(res_files)){
    if(j==1){
      system(paste("less",res_files[j],">",res_file))
    }
    if(j>1){
      system(paste("less",res_files[j],"| grep -v SNP >>",res_file))
    }
  }
}
for(num_pcs in 0:5){
  res_files = paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chrs,".assoc.logistic",sep="")
  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
  res_files = paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chrs,".assoc.logistic",sep="")
  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
  res_files = paste(curr_path,"gp_gwas_res_pcs",num_pcs,"_",chrs,".assoc.logistic",sep="")
  res_file = paste(curr_path,"gp_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
}

get_lambda_from_log_file<-function(f){
  l = readLines(f)
  s = l[grepl(" = ",l)][1]
  if(length(s)==0 || is.na(s) || nchar(s)==0){return(NA)}
  s = strsplit(s,split=" = ")[[1]][2]
  s = substr(s,start = 1,stop = 5)
  return(as.numeric(s))
}
all_log_files = list.files(curr_path)
all_log_files = all_log_files[grepl("log$",all_log_files)]
all_log_files = all_log_files[grepl("pcs",all_log_files)]
setwd(curr_path)
all_lambdas = sapply(all_log_files,get_lambda_from_log_file)
for(cc in c("elite","cooper","gp")){
  for(num_pcs in 0:5){
    curr_files = all_log_files[grepl(num_pcs,all_log_files) & grepl(cc,all_log_files)]
    print(paste(cc,num_pcs,mean(all_lambdas[curr_files],na.rm=T)))
  }
}

# reformat and write to fuma files
file2pvals = c()
for(num_pcs in 0:5){

  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F)
  res = res[,c("SNP","P")]
  colnames(res) = c("rsID","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  p = res[,2];names(p) = res[,1]
  file2pvals[[res_file]] = p

  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F)
  res = res[,c("SNP","P")]
  colnames(res) = c("rsID","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  p = res[,2];names(p) = res[,1]
  file2pvals[[res_file]] = p

  res_file = paste(curr_path,"gp_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_gp_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F)
  res = res[,c("SNP","P")]
  colnames(res) = c("rsID","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  p = res[,2];names(p) = res[,1]
  file2pvals[[res_file]] = p
  
}

table_after_name_matching <-function(x1,x2,thr){
  inds = intersect(names(x1),names(x2))
  table(x1[inds] < thr, x2[inds]<thr)
}
table_after_name_matching(file2pvals[[9]],file2pvals[[7]],1e-6)

m = c()
for(cc in c("elite","cooper","gp")){
  for(num_pcs in 0:5){
    curr_file1 = all_log_files[grepl(paste("pcs",num_pcs,sep=""),all_log_files) & grepl(cc,all_log_files)]
    curr_file2 = names(file2pvals)[grepl(paste("pcs",num_pcs,sep=""),names(file2pvals)) & grepl(cc,names(file2pvals))]
    currp = file2pvals[[curr_file2]]
    m = rbind(m,c(cc,num_pcs,all_lambdas[curr_file1][1],
          median(all_lambdas[curr_file1],na.rm=T),sum(currp < 1e-8,na.rm = T),
          median(currp,na.rm=T)))
  }
}
write.table(m,quote=F,row.names=F,sep="\t")

rownames(res_cooper) = res_cooper$SNP
rownames(res_elite) = res_elite$SNP
rownames(res_gp) = res_gp$SNP
snps = intersect(res_gp$SNP,res_elite$SNP)
snps = intersect(snps,res_cooper$SNP)
x1 = res_cooper[snps,"P"]
x2 = res_elite[snps,"P"]
x3 = res_gp[snps,"P"]
table(x1<1e-5,x2<1e-5)
table(x3<1e-5,x2<1e-5)

d = read.table(paste(curr_path,"our_subjects_covs.phe",sep=""),header=T)
pc_rocs = c()
for(j in 1:40){
  curr_inds = !is.na(d[,"elite_col"]) 
  p1 = compute_pc_vs_binary_variable_association_roc(
    pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"elite_col"]
  )
  curr_inds = !is.na(d[,"cooper_col"]) 
  p2 = compute_pc_vs_binary_variable_association_roc(
    pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"cooper_col"]
  )
  pc_rocs = rbind(pc_ps,c(p1,p2))
}
pc_qs = apply(pc_ps,2,p.adjust)
pc_inds = pc_qs < 0.01

#####################################################################################
#####################################################################################
#####################################################################################

# Take a specific PC after manual inspection and run a GWAS against it
PC = "PC6"
curr_path = paste(out_path,"gwas_",PC,"/",sep="")
system(paste("mkdir",curr_path))
pca_data = read.table(paste(out_path,"run_pca.eigenvec",sep=''))
colnames(pca_data) = c("FID","IID",paste("PC",1:40,sep=""))
write.table(pca_data,paste(curr_path,"pcs.phe",sep=""),
            row.names = F,col.names = T,quote=F,sep="\t")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
  "--linear hide-covar",
  paste("--pheno",paste(curr_path,"pcs.phe",sep="")),
  "--pheno-name",PC,
  "--threads",8,
  "--out",paste(curr_path,"gwas_res",sep="")
)
run_plink_command(curr_cmd,curr_path,"pc_gwas_run",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
wait_for_job(120)
# read the results, look at the selected SNPs
res = read.table(paste(curr_path,"gwas_res.assoc.linear",sep=""),
                 stringsAsFactors = F,row.names=2,header=T)
curr_PC_snps = rownames(res)[res$P < 1e-8]
for(chr in chrs){
  f1 = paste(out_path,"dataset1_",chr,".frq",sep='')
  f2 = paste(out_path,"dataset2_",chr,".frq",sep='')
  freqs1 = read.table(f1,header=T)
  freqs2 = read.table(f2,header=T)
  rownames(freqs1) = freqs1$SNP
  rownames(freqs2) = freqs2$SNP
  inds = intersect(freqs1$SNP,freqs2$SNP)
  inds = intersect(inds,curr_PC_snps)
  break
}

# are these directly genotypes?
direct_g = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data.bim",
                      row.names = 2,stringsAsFactors = F)
intersect(curr_PC_snps,rownames(direct_g))


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Locally
# 1000G direct analysis
setwd("/Users/David/Desktop/elite/november2018_analysis/qc_1000g_mega_direct/")
d2 = read.table(
  "/Users/David/Desktop/elite/november2018_analysis/qc_1000g/integrated_call_samples_v3.20130502.ALL.panel"
  ,stringsAsFactors=F,header=T)

# UKBB direct analysis
setwd("/Users/David/Desktop/elite/november2018_analysis/qc_ukbb_direct/")
d2 = read.table("20k_rand_controls_sex_age_with_info.txt",stringsAsFactors=F,header=F,row.names = 1)

source("~/Desktop/repos/fitness_genetics/R/gwas_flow_helper_functions.R")

pcax = read_pca_res("run_pca.eigenvec")
# d = read.table("/Users/David/Desktop/elite/november2018_analysis/qc_1000g/all_cohorts.phe",header=T,stringsAsFactors = F)
d = read.table("../qc_1000g_mega_direct/integrated_sample_metadata_and_covariates.phe",header=T,stringsAsFactors = F)
# d = read.table("../qc_1000g_mega_direct/our_subjects_covs.phe",header=T,stringsAsFactors = F)
rownames(d) = d$IID
inds = intersect(rownames(d),rownames(pcax))
d = d[inds,]
dim(d)

CohortName = d$Cohort
table(CohortName)
CohortName[CohortName=="1"]="Cooper"
CohortName[CohortName=="2"]="ELITE"
d = cbind(d,CohortName)
d$CohortName = as.character(d$CohortName)

# Other cohorts: 1000G
cohorts1 = c(d$CohortName,d2[,2])
names(cohorts1) = c(rownames(d),d2[,1])
cohorts2 = c(d$CohortName,d2[,3])
names(cohorts2) = c(rownames(d),d2[,1])
cohorts3 = c(d$CohortName,rep("1000g",nrow(d2)))
names(cohorts3) = c(rownames(d),d2[,1])
pcax = pcax[names(cohorts1),]

# Cohorts for analysis using our_subjects_covs.phe
cohorts_new = d[,c("cooper_col","elite_col","gp_col")]
cohorts_new[is.na(cohorts_new)] = "NA"
cohorts_new[cohorts_new[,1]==1,1] = "cooper"
cohorts_new[cohorts_new[,1]==2,1] = "1000g"
cohorts_new[cohorts_new[,2]==1,1] = "elite"
cohorts_new[cohorts_new[,3]==1,1] = "gp"
cohorts_new = cohorts_new[,1]; table(cohorts_new)
cohorts_new = c(cohorts_new,d2[,2])
cohorts_new [cohorts_new=="NA"]="gp"
names(cohorts_new) = c(rownames(d),d2[,1])
pcax = d[,paste("PC",1:40,sep="")]

# Other cohorts: UKBB
cohorts_all = d$CohortName
names(cohorts_all) = rownames(d)
cohorts_all[rownames(d2)]="UKBB"
table(cohorts_all)
length(intersect(rownames(pcax),names(cohorts_all)))

# PCA plots
inds = rownames(pcax)
curr_cohorts = cohorts_new
curr_cohorts = cohorts3
curr_cohorts = cohorts_all

# inds = names(cohorts2)[cohorts2=="elite" | cohorts2=="cooper" | cohorts2=="EUR"]
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],curr_cohorts[inds],curr_cohorts[inds],
      main = "All filters",xlab="PC1",ylab="PC2",lwd=3,cex=0.1)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,3],pcax[inds,4],curr_cohorts[inds],curr_cohorts[inds],
      main = "All filters",xlab="PC3",ylab="PC4",lwd=2,cex=0.1)
legend(x="bottomright",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,5],pcax[inds,6],curr_cohorts[inds],curr_cohorts[inds],
      main = "All filters",xlab="PC5",ylab="PC6",lwd=3,cex=0.01)
legend(x="topright",names(res[[1]]),fill = res[[1]])
table(pcax[,5] > 0.05)

res = two_d_plot_visualize_covariate(pcax[inds,7],pcax[inds,8],curr_cohorts[inds],curr_cohorts[inds],
    main = "All filters",xlab="PC7",ylab="PC8",lwd=2,cex=0.1)
legend(x="topright",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,9],pcax[inds,10],curr_cohorts[inds],curr_cohorts[inds],
    main = "All filters",xlab="PC9",ylab="PC10",lwd=2,cex=0.1)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,11],pcax[inds,12],curr_cohorts[inds],curr_cohorts[inds],
    main = "All filters",xlab="PC11",ylab="PC12",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,13],pcax[inds,14],curr_cohorts[inds],curr_cohorts[inds],
    main = "All filters",xlab="PC13",ylab="PC14",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# Define EUR-descendants manually
yy = pcax[,1] < 0.005 & pcax[,2]< 0.03 & pcax[,2] > -0.03
table(yy,curr_cohorts[names(yy)])

res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],yy[inds],yy[inds],
  main = "All filters",xlab="PC1",ylab="PC2",lwd=2,cex=0.1)
legend(x="topright",names(res[[1]]),fill = res[[1]])


# Do the pcs perfectly predict the group?
y = as.factor(curr_cohorts=="ukbb")
y = as.factor(cohorts3=="elite")
y = as.factor(cohorts3=="cooper")
y = as.factor(cohorts3=="1000g")
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

res_cooper = read.table("fuma_cooper_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)
res_elite = read.table("fuma_elite_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)
res_gp = read.table("fuma_gp_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)

inds = intersect(rownames(res_cooper),rownames(res_elite))
qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_elite$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)

inds = intersect(inds,rownames(res_gp))
qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_elite$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_gp$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)

plot(y=-log(res_gp[inds,"P.value"],10),x=-log(res_elite[inds,"P.value"],10),pch=20,
       ylab="GP p-value",xlab="ELITE p-value");abline(0,1)

table(res_gp[inds,"P.value"] < 10^-5,res_elite[inds,"P.value"]<10^-5 )
selected_snps = inds[res_gp[inds,"P.value"] < 10^-5]
selected_snps = selected_snps[!is.na(selected_snps)]
write.table(t(t(selected_snps)),file="gp_significant_snps_vs_1000g.txt",row.names = F,col.names = F,quote = F)

# qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
#        ylab="Sample quantiles",xlab="Theoretic quantiles",
#        xlim=c(0,20),ylim=c(0,20));abline(0,1)

# library(lattice);
# x=res$P
# qqmath(~-log10(na.omit(x)),
#        distribution=function(x){-log10(qunif(1-x))}
# )
# abline(0,1,add=T)









































