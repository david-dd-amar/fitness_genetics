
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

bfiles = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/with_ukbb/"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/integrated_sample_metadata_and_covariates.phe"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"
our_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
try({system(paste("mkdir",out_path))})
cohorts_to_exclude = ""
maf_threshold = 0.05

# The steps of the analysis below

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
covars = cbind(our_covars[,c("FID","IID","sex","batch","age","Cohort")],cohorts)
colnames(covars) = colnames(external_covars)
covars = as.matrix(rbind(covars,external_covars))
rownames(covars) = covars[,"IID"]
covars[covars[,7]=="2",7] = "elite"
covars[covars[,7]=="1",7] = "cooper"
for(j in 1:ncol(covars)){
  covars[,j] = gsub(" ","",as.character(covars[,j]))
}
covars = covars[!is.element(covars[,"CohortName"],set=cohorts_to_exclude),]

ukbb_or_not = covars[,"CohortName"] == "ukbb"
ukbb_or_not = as.numeric(ukbb_or_not)+1
covars = cbind(covars,ukbb_or_not)
colnames(covars)

write.table(file=paste(out_path,"all_cohorts.phe",sep=''),
            covars,sep=" ",row.names = F,col.names = T,quote=F)

sample_file = paste(out_path,"analysis_samples.txt",sep='')
write.table(file=sample_file,covars[,1:2],sep="\t",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# Preprocessing
# Flipscan each chromosome
for(j in 1:22){
  curr_name = paste("chr",j,"_flipscan",sep="")
  curr_cmd = paste("plink --bfile",paste(bfiles,"chr",j,sep=''),
                   "--flip-scan --allow-no-sex",
                   "--pheno",paste(out_path,"all_cohorts.phe",sep=''),
                   "--pheno-name ukbb_or_not",
                   "--keep",sample_file,
                   "--threads 4",
                   "--out",paste(out_path,curr_name,sep="")
  )
  run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(waittime = 120)
# Read the flipscan results and get the problematic SNPs
flipscan_snps = c()
for(j in 1:22){
  curr_name = paste(out_path,"chr",j,"_flipscan.flipscan",sep="")
  curr_tmp = paste(out_path,"tmp",j,sep="")
  system(paste("less",curr_name, "| grep -v -P 'NA\\s+$' >",curr_tmp))
  curr_flipscan_res = read.table(curr_tmp,stringsAsFactors = F,header=T)
  print(dim(curr_flipscan_res))
  flipscan_snps = c(flipscan_snps,curr_flipscan_res$SNP)
}
# print results to file
write.table(t(t(flipscan_snps)),file=paste(out_path,"flipscan_snps.txt",sep=""),
            row.names = F,quote = F,col.names = F)
# clean dir
system(paste("rm ",out_path,"*_flipscan*",sep=""))
system(paste("rm ",out_path,"*tmp*",sep=""))

# create a sample file for each cohort
cohorts = covars[,"CohortName"]
cohort2sample_file = c()
for(cc in unique(cohorts)){
  currfam = covars[cohorts==cc,1:2]
  curr_sample_file = paste(out_path,cc,"_samples.txt",sep="")
  cohort2sample_file[cc] = curr_sample_file
  write.table(currfam,curr_sample_file,sep="\t",quote = F,row.names = F,col.names = F)
}
for(cc in unique(cohorts)){
  curr_sample_file = paste(cc,"_samples.txt",sep="")
  for(j in 1:22){
    curr_name = paste(cc,"_chr",j,sep="")
    curr_cmd = paste("plink --bfile",paste(bfiles,"chr",j,sep=''),
                     "--keep",paste(out_path,curr_sample_file,sep=""),
                     "--maf", maf_threshold,
                     "--freq --threads 4",
                     "--out",paste(out_path,curr_name,sep="")
    )
    run_plink_command(curr_cmd,out_path,curr_name,
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
  }
}
wait_for_job(waittime = 120)
# Read all resulting freq files and clean the dir
cohort2snps = list()
for(cc in unique(cohorts)){
  curr_snps = c()
  for(j in 1:22){
    print(j)
    curr_name = paste(out_path,cc,"_chr",j,sep="")
    curr_freqs = read.table(paste(curr_name,".frq",sep=""),header=T,stringsAsFactors = F)
    curr_freqs = curr_freqs[curr_freqs$MAF > maf_threshold,]
    curr_snps = c(curr_snps,curr_freqs$SNP)
    # system(paste("rm ",curr_name,"*",sep=""))
  }
  cohort2snps[[cc]] = curr_snps
}
sapply(cohort2snps,length)
intr = cohort2snps[[1]]
for(j in 2:length(cohort2snps)){
  intr = intersect(intr,cohort2snps[[j]])
}
write.table(t(t(intr)),file=paste(out_path,"maf_filter_snps_",maf_threshold,".txt",sep=""),
            row.names = F,quote = F,col.names = F)
# clean the dir
# Read all resulting freq files and clean the dir
for(cc in unique(cohorts)){
  for(j in 1:22){
    curr_name = paste(out_path,cc,"_chr",j,sep="")
    system(paste("rm ",curr_name,"*",sep=""))
  }
}

flipscan_snp_file = paste(out_path,"flipscan_snps.txt",sep="")
maf_snp_file = paste(out_path,"maf_filter_snps_",maf_threshold,".txt",sep="")

# Add another comparison: look at the MAFs in the controls, take only
# variants with less than 1% diff in their MAF
# Look at the freqplot data - comparing mafs to 1000G
ukbb_1000g_maf_comp = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/FreqPlot-merged_control_geno-1000G.txt",
                                 stringsAsFactors = F)
diffs = abs(ukbb_1000g_maf_comp[,4])
rownames(ukbb_1000g_maf_comp) = ukbb_1000g_maf_comp[,1]
table(diffs > 0.01)
quantile(abs(ukbb_1000g_maf_comp[res_sig_snps,4]))
low_diff_snps = ukbb_1000g_maf_comp[diffs < 0.01,1]
intr = read.table(maf_snp_file,stringsAsFactors = F)[,1]
new_intr = intersect(low_diff_snps,intr)
write.table(t(t(new_intr)),file=maf_snp_file,row.names = F,quote = F,col.names = F)

####################################################################################################
####################################################################################################
####################################################################################################
# Data filtering and analysis
flipscan_snp_file = paste(out_path,"flipscan_snps.txt",sep="")
maf_snp_file = paste(out_path,"maf_filter_snps_",maf_threshold,".txt",sep="")

# LD-prune each chromosome
for(j in 1:22){
  curr_name = paste("chr",j,"_ld_prune",sep="")
  curr_cmd = paste("plink --bfile",paste(bfiles,"chr",j,sep=''),
                   "--indep-pairwise 500 10",0.1,
                   "--keep",sample_file,
                   "--exclude",flipscan_snp_file,
                   "--extract",maf_snp_file,
                   "--maf", maf_threshold,
                   "--threads 4",
                   "--out",paste(out_path,curr_name,sep="")
  )
  run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(waittime = 120)

# Sanity check: overlap between flipscan and pruned list?
x1 = read.table(flipscan_snp_file,stringsAsFactors = F)[,1]
x2 = read.table(paste(out_path,"chr2_ld_prune_filter.bim",sep=""),stringsAsFactors = F)[,2]
length(intersect(x1,x2))

# For each chromosome keep ld pruned snps and save bed
out_pruned_beds = c()
for(j in 1:22){
  curr_name = paste("chr",j,"_ld_prune_filter",sep="")
  curr_snps = paste(out_path,"chr",j,"_ld_prune.prune.in",sep="")
  curr_cmd = paste("plink --bfile",paste(bfiles,"chr",j,sep=''),
                   "--keep",sample_file,
                   "--extract", curr_snps,
                   "--threads 4",
                   "--make-bed --out",paste(out_path,curr_name,sep="")
  )
  out_pruned_beds[j] = paste(out_path,curr_name,sep="")
  run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(waittime = 120)
# Clean the directory
setwd(out_path)
system(paste("rm","*.out"))
system(paste("rm","*.nosex"))
system(paste("rm","*.err"))
system(paste("rm","*.flipscan"))

# Merge the beds before PCA and relatedness analyses
allfiles_path = paste(out_path,"out_pruned_beds.txt",sep="")
write.table(t(t(out_pruned_beds[-1])),
            file = allfiles_path,sep="",row.names = F,col.names = F,quote = F)
curr_cmd = paste("plink --bfile",out_pruned_beds[1],
                 "--merge-list",allfiles_path,
                 "--threads 4",
                 "--make-bed --out",paste(out_path,"merged_ld_pruned",sep='')
)
run_plink_command(curr_cmd,out_path,"merged_ld_pruned",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)

####################################################################################################
####################################################################################################
####################################################################################################
# Population stratification analysis
# We test several methods, these will hopefully help in understanding the intra-EU
# structure that is the best one for doing the GWAS.
# Tested methods:
# PLINK PCA, PCAngsd admix, fastStructure, ADMIXTURE
# Each method is run before and after the manual selection of EU samples
# we put all data in the same dir
pop_anal_dir = paste(out_path,"pop_analysis/",sep="")
system(paste("mkdir",pop_anal_dir))
# Load clustering of the data (manual or by using a differnt method/script)
load(paste(out_path,"eu_clustering_elite_data_manual.RData",sep=""))
bfile_all = paste(out_path,"merged_ld_pruned",sep='')
bfile_eu = paste(out_path,"filters_cleaned_data",sep='')
bfiles_pop = list("all" = bfile_all,"eu"=bfile_eu)

curr_cmd = paste("plink --bfile",bfile_all,
                 "--threads 4",
                 "--keep",paste(out_path,"filters_cleaned_subjects.txt",sep=''),
                 "--make-bed",
                 "--out",paste(out_path,"filters_cleaned_data",sep='')
)
run_plink_command(curr_cmd,out_path,"filtered_bysample_dataset",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
wait_for_job(waittime = 120)

for(nn in names(bfiles_pop)){
  # Run PCA and relatedness 
  # In theory we should remove rl subjects but we assume that this is not a major issue.
  # In addition, standard filters will remove most Jap subjects which we do not want for
  # a global analysis.
  curr_name = paste("pca_rl_",nn,sep='')
  curr_cmd = paste("plink --bfile",bfiles_pop[[nn]],
                   "--genome --min 0.2",
                   "--threads 16",
                   "--pca 40",
                   "--out",paste(pop_anal_dir,curr_name,sep='')
  )
  run_plink_command(curr_cmd,pop_anal_dir,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)
  wait_for_job(waittime = 120)
  
  # PCANgsd
  setwd(pop_anal_dir)
  err_path = paste(pop_anal_dir,"pcangsd_run_",nn,".err",sep="")
  log_path = paste(pop_anal_dir,"pcangsd_run_",nn,".log",sep="")
  curr_cmd = paste("module load py-scipystack\n",
                   "python /home/users/davidama/repos/pcangsd/pcangsd.py",
                   "-plink",bfile_all,
                   "-admix",
                   "-threads 16",
                   "-o",paste(pop_anal_dir,"pcangsd_",nn,sep=''))
  curr_sh_file = paste(pop_anal_dir,"pcangsd_run_",nn,".sh",sep="")
  print_sh_file(curr_sh_file,
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=16,mem_size=64000,time="12:00:00"),curr_cmd)
  system(paste("sbatch",curr_sh_file))
  wait_for_job(waittime = 120)
  
  # fastStructure
  faststruct_name = paste("faststructure_",nn,sep="")
  faststruct_path = paste(pop_anal_dir,faststruct_name,"/",sep="")
  system(paste("mkdir",faststruct_path))
  setwd(faststruct_path)
  for (k in 1:10){
    # 1. Run the algorithm on a number of Ks
    err_path = paste("run_",k,".err",sep="")
    log_path = paste("run_",k,".log",sep="")
    curr_cmd = paste("module load py-scipystack\n",
                     "module load py-cython\n",
                     "module load gsl\n",
                     "python /home/users/davidama/apps/fastStructure/structure.py",
                     "-K",k,
                     paste("--input=",bfiles_pop[[nn]],sep=''),
                     "--output=model",
                     "--full --seed=123 --tol=10e-5"
    )
    curr_sh_file = paste("run_",k,".sh",sep="")
    print_sh_file(curr_sh_file,
                  get_sh_prefix_one_node_specify_cpu_and_mem(
                    err_path,log_path,Ncpu=4,mem_size=32000,time="12:00:00"),curr_cmd)
    system(paste("sbatch",curr_sh_file))  
  }
  
  # ADMIXTURE
  adm_name = paste("admixture_",nn,sep="")
  adm_path = paste(pop_anal_dir,adm_name,"/",sep="")
  system(paste("mkdir",adm_path))
  setwd(adm_path)
  for (k in 1:10){
    # 1. Run the algorithm on a number of Ks
    err_path = paste("run_",k,".err",sep="")
    log_path = paste("run_",k,".log",sep="")
    curr_cmd = paste(
                     "/home/users/davidama/apps/admixture/admixture_linux-1.3.0/admixture",
                     paste(bfiles_pop[[nn]],".bed",sep=""), k,"--cv",
                     "-s 123 -j8 -C 0.1"
          )
    curr_sh_file = paste("run_",k,".sh",sep="")
    print_sh_file(curr_sh_file,
                  get_sh_prefix_one_node_specify_cpu_and_mem(
                    err_path,log_path,Ncpu=8,mem_size=32000,time="18:00:00"),curr_cmd)
    system(paste("sbatch",curr_sh_file))  
  }
}

subjects_for_analysis = rownames(d)[all_filters]
newd = cbind(d[subjects_for_analysis,],pcax[subjects_for_analysis,])
save(subjects_for_analysis,pcax,file=paste(out_path,"clustering_data.RData",sep=""))
write.table(file=paste(out_path,"covariates_filters_cleaned.phe",sep=''),
            newd,sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=paste(out_path,"filters_cleaned_subjects.txt",sep=''),
            newd[,1:2],sep="\t",row.names = F,col.names = T,quote=F)

# Write the pheno files for each analysis
subjects_for_analysis = intersect(rownames(new_pcax),rownames(d))
newd = cbind(d[subjects_for_analysis,],new_pcax[subjects_for_analysis,])
write.table(file=paste(out_path,"filters_cleaned_all.phe",sep=''),
            newd,sep=" ",row.names = F,col.names = T,quote=F)
newd1 = newd[newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="cooper",]
newd1$ExerciseGroup = newd1$ExerciseGroup-1
table(newd1$ExerciseGroup)
write.table(file=paste(out_path,"filters_cleaned_elite_vs_cooper.phe",sep=''),
            newd1,sep=" ",row.names = F,col.names = T,quote=F)
newd2 = newd[newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="ukbb",]
newd2$ExerciseGroup[newd2$ExerciseGroup > 1] = 2
table(newd2$ExerciseGroup)
write.table(file=paste(out_path,"filters_cleaned_elite_vs_ukbb.phe",sep=''),
            newd2,sep=" ",row.names = F,col.names = T,quote=F)
newd3 = newd[newd[,"CohortName"] == "cooper" | newd[,"CohortName"]=="ukbb",]
newd3$ExerciseGroup[newd3$ExerciseGroup > 1] = 2
table(newd3$ExerciseGroup)
write.table(file=paste(out_path,"filters_cleaned_cooper_vs_ukbb.phe",sep=''),
            newd3,sep=" ",row.names = F,col.names = T,quote=F)

# # Write the pheno files for each analysis
# new_pcax = read.table(paste(out_path,"merged_ld_pruned_pcangsd.K16.a0.qopt",sep=""))
# fam = read.table(paste(out_path,"filters_cleaned_data.fam",sep=""),stringsAsFactors = F)
# rownames(new_pcax) = fam[,2]
# new_pcax = new_pcax[,-ncol(new_pcax)] # last col is linearly dependent on the prev ones
# colnames(new_pcax) = paste("PC",1:ncol(new_pcax),sep="") # make it in the expected format below
# d = read.table(paste(out_path,"all_cohorts.phe",sep=''),header=T,stringsAsFactors = F)
# rownames(d) = d$IID
# subjects_for_analysis = intersect(rownames(new_pcax),rownames(d))
# newd = cbind(d[subjects_for_analysis,],new_pcax[subjects_for_analysis,])
# 
# write.table(file=paste(out_path,"pcangsd_filters_cleaned_all.phe",sep=''),
#             newd,sep=" ",row.names = F,col.names = T,quote=F)
# newd1 = newd[newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="cooper",]
# newd1$ExerciseGroup = newd1$ExerciseGroup-1
# table(newd1$ExerciseGroup)
# write.table(file=paste(out_path,"pcangsd_filters_cleaned_elite_vs_cooper.phe",sep=''),
#             newd1,sep=" ",row.names = F,col.names = T,quote=F)
# newd2 = newd[newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="ukbb",]
# newd2$ExerciseGroup[newd2$ExerciseGroup > 1] = 2
# table(newd2$ExerciseGroup)
# write.table(file=paste(out_path,"pcangsd_filters_cleaned_elite_vs_ukbb.phe",sep=''),
#             newd2,sep=" ",row.names = F,col.names = T,quote=F)
# newd3 = newd[newd[,"CohortName"] == "cooper" | newd[,"CohortName"]=="ukbb",]
# newd3$ExerciseGroup[newd3$ExerciseGroup > 1] = 2
# table(newd3$ExerciseGroup)
# write.table(file=paste(out_path,"pcangsd_filters_cleaned_cooper_vs_ukbb.phe",sep=''),
#             newd3,sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# GWAS

flipscan_snp_file = paste(out_path,"flipscan_snps.txt",sep="")
maf_snp_file = paste(out_path,"maf_filter_snps_",maf_threshold,".txt",sep="")

# Run all GWAS tests
covar_files = c(
  paste(out_path,"filters_cleaned_elite_vs_cooper.phe",sep=''),
  paste(out_path,"filters_cleaned_elite_vs_ukbb.phe",sep=''),
  paste(out_path,"filters_cleaned_cooper_vs_ukbb.phe",sep=''),
  paste(out_path,"pcangsd_filters_cleaned_elite_vs_cooper.phe",sep=''),
  paste(out_path,"pcangsd_filters_cleaned_elite_vs_ukbb.phe",sep=''),
  paste(out_path,"pcangsd_filters_cleaned_cooper_vs_ukbb.phe",sep=''),
  paste(out_path,"filters_cleaned_elite_vs_cooper.phe",sep=''),
  paste(out_path,"filters_cleaned_elite_vs_ukbb.phe",sep=''),
  paste(out_path,"filters_cleaned_cooper_vs_ukbb.phe",sep='')
)
PCs = list(
  paste("PC",1:6,sep=""),
  paste("PC",1:8,sep=""),
  paste("PC",c(1:9,14),sep=""),
  paste("PC",1:15,sep=""),
  paste("PC",1:15,sep=""),
  paste("PC",1:15,sep=""),
  "",
  "",
  ""
)
gwas_dir = paste(out_path,"gwas/",sep="")
system(paste("mkdir",gwas_dir))

analysis_names = c(rep(NA,6),
                   "filters_cleaned_elite_vs_cooper_0PCs",
                   "filters_cleaned_elite_vs_ukbb_0PCs",
                   "filters_cleaned_cooper_vs_ukbb_0PCs")

for(i in 7:length(covar_files)){
  covar_file = covar_files[i]
  curr_name = analysis_names[i]
  if(is.na(curr_name)){
    curr_name = gsub(covar_file,pattern = ".phe",replacement = "")
    curr_name = strsplit(curr_name,split="/")[[1]]
    curr_name = curr_name[length(curr_name)] 
  }
  curr_dir = paste(gwas_dir,curr_name,"/",sep="")
  system(paste("mkdir",curr_dir))
  curr_PCs = PCs[[i]]
  for(j in 1:22){
    covars_line = paste("--covar-name sex,Age,",paste(curr_PCs,collapse=","),sep="")
    if(length(curr_PCs)==1 && curr_PCs==""){
      covars_line = "--covar-name sex,Age"
    }
    curr_cmd = paste("plink2",
                     "--bfile",paste(bfiles,"chr",j,sep=''),
                     "--logistic hide-covar firth-fallback",
                     "--exclude",flipscan_snp_file,
                     "--extract",maf_snp_file,
                     "--threads 4",
                     paste("--pheno",covar_file),
                     paste("--pheno-name ExerciseGroup"),
                     paste("--covar",covar_file),covars_line,
                     "--adjust --out",paste(curr_dir,"chr",j,sep=''))
    run_plink_command(curr_cmd,curr_dir,paste("logistic_chr",j,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000,
                      plink_pkg = "plink/2.0a1")
  }
}

# Read some results
res = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/gwas/filters_cleaned_cooper_vs_ukbb_4PCs/chr1.ExerciseGroup.glm.logistic.hybrid.adjusted")
quantile(res[,3])

dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/gwas/filters_cleaned_cooper_vs_ukbb_0PCs/"
curr_adj_files = list.files(dir)
curr_adj_files = curr_adj_files[grepl("adjusted$",curr_adj_files)]
for(ff in curr_adj_files){
  print(ff)
  res = read.table(paste(dir,ff,sep=""),stringsAsFactors = F)
  print(res[1:2,1:3])
}

dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/gwas/filters_cleaned_elite_vs_ukbb_0PCs/"
curr_log_files = list.files(dir)
curr_log_files = curr_log_files[grepl("log$",curr_log_files) & grepl("logistic",curr_log_files)]
for(ff in curr_log_files){
  print(ff)
  res = readLines(paste(dir,ff,sep=""))
  nres = length(res)
  print(res[(nres-5):nres])
}


# Add some PC-GWASs
all_covars_file = paste(out_path,"filters_cleaned_all.phe",sep='')
for(i in c(2:3,6)){
  covar_file = all_covars_file
  curr_name = paste("PC",i,sep="")
  curr_dir = paste(gwas_dir,curr_name,"/",sep="")
  system(paste("mkdir",curr_dir))
  for(j in 1:22){
    curr_cmd = paste("plink2",
                     "--bfile",paste(bfiles,"chr",j,sep=''),
                     "--linear hide-covar",
                     "--exclude",flipscan_snp_file,
                     "--extract",maf_snp_file,
                     "--threads 4",
                     "--pheno",covar_file,
                     "--pheno-name", curr_name,
                     paste("--covar",covar_file),
                     "--covar-name sex,Age",
                     "--adjust --out",paste(curr_dir,"chr",j,sep=''))
    run_plink_command(curr_cmd,curr_dir,paste("linear_chr",j,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000,
                      plink_pkg = "plink/2.0a1")
  }
}
wait_for_job(waittime = 120)

####################################################################################################
####################################################################################################
####################################################################################################
# QA: check if JHU SNPs separate the groups
load("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/with_ukbb/bim_overlap_analysis_results.RData")
table(grepl("JHU",final_shared_snps[,1]))
jhu_snps = final_shared_snps[grepl("JHU",final_shared_snps[,2]),1]

dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/gwas/PC6/"
curr_adj_files = list.files(dir)
curr_adj_files = curr_adj_files[grepl("adj",curr_adj_files)]
for(ff in curr_adj_files){
  print(ff)
  res = read.table(paste(dir,ff,sep=""),stringsAsFactors = F,header=T)
  rownames(res) = res[,1]
  is_jhu = is.element(res[,1],set=jhu_snps)
  x1 = res[is_jhu,2]
  x2 = res[!is_jhu,2]
}

dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/gwas/filters_cleaned_cooper_vs_ukbb_4PCs/"
curr_adj_files = list.files(dir)
curr_adj_files = curr_adj_files[grepl("adjusted$",curr_adj_files)]
for(ff in curr_adj_files){
  res = read.table(paste(dir,ff,sep=""),stringsAsFactors = F)
  rownames(res) = res[,2]
  is_jhu = is.element(res[,2],set=jhu_snps)
  x1 = res[is_jhu,3]
  x2 = res[!is_jhu,3]
  break
}




####################################################################################################
####################################################################################################
####################################################################################################
# # Look at the freqplot data - comparing mafs to 1000G
# ukbb_1000g_maf_comp = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/FreqPlot-merged_control_geno-1000G.txt",
#                                  stringsAsFactors = F)
# diffs = abs(ukbb_1000g_maf_comp[,4])
# rownames(ukbb_1000g_maf_comp) = ukbb_1000g_maf_comp[,1]
# table(diffs > 0.01)
# quantile(abs(ukbb_1000g_maf_comp[res_sig_snps,4]))
# Examples for possible strand flip: rs4970408,rs2072926,rs10909801 (all imputed)
# Manual checks of these:
# rs10909801
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/chr1.bim | grep rs10909801
# 1	rs10909801:3759937:A:G	0	3759937	A	G
# less ../../../mega_eu_imp/impute2_1000gRef_out/chr1.frq | grep 3759937
# 1  rs10909801:3759937:A:G    G    A       0.4299     3610
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/merged_geno.bim | grep rs10909801
# 1	rs10909801:3759937:A:G	0	3759937	A	G
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/merged_geno.frq | grep rs10909801
# UKBB:
# less ../../../../ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq | grep rs10909801
# rs10909801    A    G        0.485    39900
# less ../../../../ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.bim | grep rs10909801
# 1	rs10909801	0	3759937	A	G
# Summary: our dataset shows G-A for the SNP, ukbb shows A G, which is consistent with 1000G
#
#
# rs2072926
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/chr1.bim | grep rs2072926
# 1	rs2072926:1684758:G:A	0	1684758	G	A
# less ../../../mega_eu_imp/impute2_1000gRef_out/chr1.frq | grep 1684758
# 1 :1684758:G:A    A    G       0.4957     2822
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/merged_geno.bim | grep rs2072926
# 1	rs2072926:1684758:G:A	0	1684758	G	A
# less ../../../mega_eu_imp/impute2_1000gRef_out/check_bim_res/merged_geno.frq | grep rs2072926
# UKBB:
# less ../../../../ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq | grep rs2072926
# rs2072926    A    G       0.4933    39884
# less ../../../../ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.bim | grep rs2072926
# 1	rs2072926	0	1684758	G	A
# Summary: our dataset shows A G for the SNP, ukbb shows A G, which is consistent with 1000G
# IMPORTANT: how things look like before the merge in plink:
# less ../mega_eu_imp/with_ukbb/new_bed_1.frq | grep rs2072926
# rs2072926    G    A       0.5067
# less ../mega_eu_imp/with_ukbb/new_bed_2.frq | grep rs2072926
# rs2072926:1684758:G:A    G    A       0.5043
#
# Direct analysis of raw genotypes
# setwd("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_vs_ukbb_20k/qa_rs2072926")
# x0 = read.table("tmp.ped")
# rownames(x0) = x0[,2]
# x1 = read.table("tmp1.ped")
# rownames(x1) = x1[,2]
# x2 = read.table("tmp2.ped")
# rownames(x2) = x2[,2]
# table(x0[rownames(x2),7:8]==x2[,7:8])

####################################################################################################
####################################################################################################
####################################################################################################

# Merge all results to a single file
gwas_dir = paste(out_path,"gwas/",sep="")
out_folders = list.files(gwas_dir)
for(i in 1:length(out_folders)){
  curr_dir = paste(gwas_dir,out_folders[i],"/",sep="")
  print(curr_dir)
  curr_files = list.files(curr_dir)
  setwd(curr_dir)
  
  res_files1 = curr_files[grepl("logistic.hybrid$",curr_files) |
                            grepl("glm.linear$",curr_files)]
  res_file1 = paste(curr_dir,"plink_all_assoc.txt",sep="")
  for(j in 1:length(res_files1)){
    if(j==1){
      system(paste("less",res_files1[j],">",res_file1))
    }
    if(j>1){
      system(paste("less",res_files1[j],"| grep -v CHROM >>",res_file1))
    }
  }
  
  res_files2 = curr_files[grepl("adjusted$",curr_files)]
  res_file2 = paste(curr_dir,"plink_all_assoc_adj.txt",sep="")
  for(j in 1:length(res_files2)){
    if(j==1){
      system(paste("less",res_files2[j],">",res_file2))
    }
    if(j>1){
      system(paste("less",res_files2[j],"| grep -v CHROM >>",res_file2))
    }
  }
  for(ff in union(res_files2,res_files1)){
    system(paste("rm",ff))
  }
}
for(i in 1:length(out_folders)){
  curr_dir = paste(gwas_dir,out_folders[i],"/",sep="")
  res_file2 = paste(curr_dir,"plink_all_assoc_adj.txt",sep="")
  print(res_file2)
  res = read.table(res_file2)
  # res = res[,2:3]
  colnames(res) = c("rsID","P-value")
  write.table(res,res_file2,row.names = F,col.names = T,quote=F,sep=" ")
}

####################################################################################################
####################################################################################################
####################################################################################################

# # Sanity check 2: UKBB vs. not UKBB: flip scan
# d = read.table(paste(out_path,"kmeans_cleaned.phe",sep=''),header=T,stringsAsFactors = F)
# rownames(d) = d$IID
# covars_copy = d
# covars_copy$ExerciseGroup[covars_copy$CohortName=="ukbb"]="1"
# covars_copy$ExerciseGroup[covars_copy$CohortName!="ukbb"]="2"
# table(covars_copy$ExerciseGroup)
# covar_file = paste(out_path,"ukbb_vs_nonukbb.phe",sep='')
# write.table(file=covar_file,covars_copy,sep=" ",row.names = F,col.names = T,quote=F)
# err_path = paste(out_path,"ukbb_vs_nonukbb.err",sep="")
# log_path = paste(out_path,"ukbb_vs_nonukbb.log",sep="")
# curr_cmd = paste("plink",
#                  "--bfile",gwas_bfile,"--logistic --flip-scan --allow-no-sex --test-missing",
#                  paste("--pheno",covar_file),
#                  paste("--pheno-name ExerciseGroup"),
#                  "--adjust --out",paste(out_path,"ukbb_vs_nonukbb_logistic",sep=''))
# curr_sh_file = "ukbb_vs_nonukbb_logistic.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size = 10000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

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
####################################################################################################
####################################################################################################
####################################################################################################

setwd("/Users/David/Desktop/elite/november2018_analysis/mega_vs_ukbb/")
pcax = read_pca_res("merged_ld_pruned.eigenvec")
pcax = read_pca_res("filters_cleaned_pca.eigenvec")
d = read.table("all_cohorts.phe",header=T,stringsAsFactors = F)
rownames(d) = d$IID
inds = intersect(rownames(d),rownames(pcax))
pcax = pcax[inds,]
d = d[inds,]
dim(d)

# Analyze the relatedness report
library("igraph")
genome_rl_file = "merged_ld_pruned.genome"
rl_data = read.table(genome_rl_file,header=T,stringsAsFactors = F)
rl_edges = as.matrix(rl_data[,c("IID1","IID2")])
mode(rl_edges) = "character"
rl_g = igraph::graph_from_edgelist(rl_edges,directed = F)
rl_clusters = clusters(rl_g)[[1]]
rl_subjects_to_remove = c()
for(cl in unique(rl_clusters)){
  curr_subjects = names(rl_clusters)[rl_clusters==cl]
  rl_subjects_to_remove = c(rl_subjects_to_remove,curr_subjects[-1])
}
rl_clustering = is.element(rownames(d),set=rl_subjects_to_remove)
names(rl_clustering) = rownames(d)
table(rl_clustering,d$CohortName)

outliers_to_rem = rep(F,nrow(d))
names(outliers_to_rem) = rownames(d)
for(j in 1:40){
  x = pcax[,j]
  x = abs(x-mean(x))/sd(x)
  outliers_to_rem[x>7] = T
}
tb = table(outliers_to_rem,d$CohortName)
for(j in 1:ncol(tb)){
  print(tb[2,j]/sum(tb[,j]))
}
table(outliers_to_rem,d$CohortName)

test_inds = 1:nrow(d)
fisher.test(table(outliers_to_rem[test_inds],
                  d$CohortName[test_inds]=="elite"),alternative = "g")$p.value
fisher.test(table(outliers_to_rem[test_inds],
                  d$CohortName[test_inds]=="cooper"),alternative = "g")$p.value

eu_af_clustering = pcax[,2] > -0.02  & pcax[,2] < 0.03 &  pcax[,1] < 0.03 
table(eu_af_clustering)
table(eu_af_clustering,d$CohortName)

all_filters = eu_af_clustering & !rl_clustering & !outliers_to_rem
names(all_filters) = rownames(d)
table(all_filters,d$CohortName)
table(all_filters)
save(eu_af_clustering,rl_clustering,outliers_to_rem,all_filters,file="eu_clustering_elite_data_manual.RData")

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],all_filters[inds],all_filters[inds],
      main = "All filters",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],eu_af_clustering[inds],eu_af_clustering[inds],
    main = "Non-Asian clustering",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],d[inds,]$CohortName,d[inds,]$CohortName,
    main = "By cohort",xlab="PC1",ylab="PC2",lwd=2)
legend(x="bottomright",names(res[[1]]),fill = res[[1]])

par(mfrow=c(1,2))
inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],rl_clustering[inds],rl_clustering[inds],
    main = "Relatedness filter",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]],cex=0.9)

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],outliers_to_rem[inds],outliers_to_rem[inds],
    main = "Outlier filter",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]],cex=0.9)

res = two_d_plot_visualize_covariate(pcax[inds,22],pcax[inds,22],outliers_to_rem[inds],outliers_to_rem[inds],
    main = "Outliers",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]])

simple_cl = abs(pcax[inds,21]-mean(pcax[inds,22]))/sd(pcax[inds,22]) > 6
res = two_d_plot_visualize_covariate(pcax[inds,21],pcax[inds,22],simple_cl,simple_cl,
    main = "Outliers",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]])

par(mfrow=c(2,2))
for(j in seq(9,16,by=2)){
  inds = rownames(d)
  res = two_d_plot_visualize_covariate(pcax[inds,j],pcax[inds,j+1],d[inds,]$CohortName,d[inds,]$CohortName,
      xlab=j,ylab=j+1,lwd=2)
  # legend(x="topleft",names(res[[1]]),fill = res[[1]],cex=0.8)
}
dev.off()
res = two_d_plot_visualize_covariate(pcax[inds,21],pcax[inds,1],d[inds,]$CohortName,d[inds,]$CohortName,
    main = "By cohort",xlab=2,ylab=3)


# Check the new PCs for association with group
pc_ps = c()
for(j in 1:40){
  curr_inds = d[,"CohortName"] == "elite" | d[,"CohortName"]=="cooper"
  p1 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  curr_inds = d[,"CohortName"] == "elite" | d[,"CohortName"]=="ukbb"
  p2 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  curr_inds = d[,"CohortName"] == "cooper" | d[,"CohortName"]=="ukbb"
  p3 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  pc_ps = rbind(pc_ps,c(p1,p2,p3))
}
colnames(pc_ps) = c("elite_vs_cooper","elite_vs_ukbb","cooper_vs_ukbb")
pc_qs = apply(pc_ps,2,p.adjust)
pc_inds = pc_qs < 0.01

# Check the new PCs for association with group
pc_rocs = c()
for(j in 1:40){
  curr_inds = d[,"CohortName"] == "elite" | d[,"CohortName"]=="cooper"
  p1 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  curr_inds = d[,"CohortName"] == "elite" | d[,"CohortName"]=="ukbb"
  p2 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  curr_inds = d[,"CohortName"] == "cooper" | d[,"CohortName"]=="ukbb"
  p3 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"CohortName"]
  )
  pc_rocs = rbind(pc_rocs,c(p1,p2,p3))
}
pc_rocs > 0.75

inds = d[,"CohortName"] == "elite" | d[,"CohortName"]=="cooper"
inds = rownames(d)
inds = d[,"CohortName"] == "cooper" | d[,"CohortName"]=="ukbb"
inds = d[,"CohortName"] == "cooper"
inds =  d[,"CohortName"]=="ukbb"
table(inds)
res = two_d_plot_visualize_covariate(pcax[inds,2],pcax[inds,6],d[inds,]$CohortName,d[inds,]$CohortName,
    main = "By cohort",xlab="PC2",ylab="PC6",lwd=2)
legend(x="bottomright",names(res[[1]]),fill = res[[1]])

# Plot fastStructure results
fam = read.table("filters_cleaned_data.fam",stringsAsFactors = F)
clust_res = read.table("5.5.meanQ")

# Plot PCANGSD results 
clust_res = read.table("merged_ld_pruned_pcangsd.K16.a0.qopt")
fam = read.table("filters_cleaned_data.fam",stringsAsFactors = F)


colSums(clust_res > 0.1)
cor(clust_res)
table(rowSums(clust_res > 0.1))
curr_clust = clust_res[,14] > 0.05
table(curr_clust)
inds = fam[,2]
names(curr_clust) = inds
table(curr_clust,d[inds,"CohortName"])
chisq.test(table(curr_clust,d[inds,"CohortName"]))$p.value
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],curr_clust,curr_clust,
    main = "By cohort",xlab="PC2",ylab="PC6",lwd=2)
res = two_d_plot_visualize_covariate(pcax[inds,3],pcax[inds,4],curr_clust,curr_clust,
    main = "By cohort",xlab="PC2",ylab="PC6",lwd=2)
res = two_d_plot_visualize_covariate(pcax[inds,5],pcax[inds,6],curr_clust,curr_clust,
    main = "By cohort",xlab="PC2",ylab="PC6",lwd=2)

ages = as.numeric(d[inds,"Age"])
cor.test(clust_res[!is.na(ages),5],ages[!is.na(ages)])

all(rownames(pcax)==fam[,1])
library(corrplot)
corrplot(cor(pcax[inds,1:20],clust_res))

yy = pcax[inds,5]
xx = as.matrix(clust_res[,1:15])
summary(lm(yy~xx))

pcax1 = read_pca_res("merged_ld_pruned.eigenvec")
pcax2 = read_pca_res("filters_cleaned_pca.eigenvec")
inds = intersect(rownames(pcax1),rownames(pcax2))
corrplot(cor(pcax1[inds,1:20],pcax2[inds,1:20]))
pcavals1 = read.table("merged_ld_pruned.eigenval")[,1]
pcavals2 = read.table("filters_cleaned_pca.eigenval")[,1]
par(mfrow=c(1,2))
plot(pcavals1,type="b",xlab="PC",ylab="Variance",pch=20,main="After filter")
plot(pcavals2,type="b",xlab="PC",ylab="Variance",pch=20,main="Before filter")
dev.off()
lambdas = c(1.6382,1.63768,1)
barplot(lambdas,beside = T,space = 0,ylim=c(1,1.7),
        xpd=F,ylab="Genomic inflation",xlab="PCs",names.arg = 4:6)
# # redu pca directly from pcax
# lambdas = read.table("merged_ld_pruned.eigenval",stringsAsFactors = F)[,1]
# xx = pcax
# for(j in 1:40){
#   xx[,j]=xx[,j]*sqrt(lambdas[j])
# }
# xx = xx[all_filters,]
# dim(xx)
# new_xx = princomp(xx,scores = T)$scores
# pcax = new_xx
# colnames(pcax) = paste("PC",1:40,sep="")

