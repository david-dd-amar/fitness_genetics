
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
maf_threshold = 0.01

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
wait_for_job(waittime = 60)
# Run PCA and relatedness
curr_cmd = paste("plink --bfile",paste(out_path,"merged_ld_pruned",sep=''),
                 "--genome --min 0.2",
                 "--pca 40",
                 "--threads 8",
                 "--out",paste(out_path,"merged_ld_pruned",sep='')
)
run_plink_command(curr_cmd,out_path,"merged_ld_pruned_pca",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
wait_for_job(waittime = 120)
# Use PCANgsd to get EU samples
setwd(out_path)
err_path = paste(out_path,"pcangsd_run.err",sep="")
log_path = paste(out_path,"pcangsd_run.log",sep="")
curr_cmd = paste("module load py-scipystack\n",
                 "python /home/users/davidama/repos/pcangsd/pcangsd.py",
                 "-plink",paste(out_path,"merged_ld_pruned",sep=''),
                 "-admix",
                 "-selection 1",
                 "-threads 16",
                 "-o",paste(out_path,"merged_ld_pruned_pcangsd",sep=''))
curr_sh_file = "pcangsd_run.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=16,mem_size=64000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(waittime = 120)

# Analyze the relatedness report
library("igraph",lib.loc = "~/R/packages")
genome_rl_file = paste(out_path,"merged_ld_pruned.genome",sep='')
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

# Load pca results and metadata
pca_file = paste(out_path,"merged_ld_pruned.eigenvec",sep='')
pcax = read_pca_res(pca_file)
d = read.table(paste(out_path,"all_cohorts.phe",sep=''),header=T,stringsAsFactors = F)
rownames(d) = d$IID
pcax = pcax[rownames(d),]
# Load clustering of the data (manual or by using a differnt method/script)
load(paste(out_path,"eu_clustering_elite_data_manual.RData",sep=""))

subjects_for_analysis = rownames(d)[all_filters]
newd = cbind(d[subjects_for_analysis,],pcax[subjects_for_analysis,])
save(subjects_for_analysis,pcax,file=paste(out_path,"clustering_data.RData",sep=""))
write.table(file=paste(out_path,"covariates_filters_cleaned.phe",sep=''),
            newd,sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=paste(out_path,"filters_cleaned_subjects.txt",sep=''),
            newd[,1:2],sep="\t",row.names = F,col.names = T,quote=F)

# Rerun PCA using the current set of subjects
curr_cmd = paste("plink --bfile",paste(out_path,"merged_ld_pruned",sep=''),
                 "--pca 40 --freq",
                 "--threads 16",
                 "--keep",paste(out_path,"filters_cleaned_subjects.txt",sep=''),
                 "--out",paste(out_path,"filters_cleaned_pca",sep='')
)
run_plink_command(curr_cmd,out_path,"eu_pruned_pca",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)
wait_for_job(waittime = 120)

pca_file = paste(out_path,"filters_cleaned_pca.eigenvec",sep='')
new_pcax = read_pca_res(pca_file)
new_pcax = new_pcax[rownames(newd),]

# Check the new PCs for association with group
pc_ps = c()
for(j in 1:40){
  curr_inds = newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="cooper"
  p1 = compute_pc_vs_binary_variable_association_p(
    pc = new_pcax[curr_inds,paste("PC",j,sep="")],y = newd[curr_inds,"CohortName"]
  )
  curr_inds = newd[,"CohortName"] == "elite" | newd[,"CohortName"]=="ukbb"
  p2 = compute_pc_vs_binary_variable_association_p(
    pc = new_pcax[curr_inds,paste("PC",j,sep="")],y = newd[curr_inds,"CohortName"]
  )
  curr_inds = newd[,"CohortName"] == "cooper" | newd[,"CohortName"]=="ukbb"
  p3 = compute_pc_vs_binary_variable_association_p(
    pc = new_pcax[curr_inds,paste("PC",j,sep="")],y = newd[curr_inds,"CohortName"]
  )
  pc_ps = rbind(pc_ps,c(p1,p2,p3))
}
colnames(pc_ps) = c("elite_vs_cooper","elite_vs_ukbb","cooper_vs_ukbb")
pc_qs = apply(pc_ps,2,p.adjust)
pc_inds = pc_ps < 0.001

# Write the pheno files for each analysis
newd = cbind(d[subjects_for_analysis,],new_pcax[subjects_for_analysis,])
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

####################################################################################################
####################################################################################################
####################################################################################################

# Run all GWAS tests
covar_files = c(
  paste(out_path,"filters_cleaned_elite_vs_cooper.phe",sep=''),
  paste(out_path,"filters_cleaned_elite_vs_ukbb.phe",sep=''),
  paste(out_path,"filters_cleaned_cooper_vs_ukbb.phe",sep='')
)
PCs = list(
  paste("PC",1:4,sep=""),
  paste("PC",1:4,sep=""),
  paste("PC",1:4,sep="")
)
gwas_dir = paste(out_path,"gwas/",sep="")
system(paste("mkdir",gwas_dir))

for(i in 1:length(covar_files)){
  covar_file = covar_files[i]
  curr_name = gsub(covar_file,pattern = ".phe",replacement = "")
  curr_name = strsplit(curr_name,split="/")[[1]]
  curr_name = curr_name[length(curr_name)]
  curr_dir = paste(gwas_dir,curr_name,"/",sep="")
  system(paste("mkdir",curr_dir))
  curr_PCs = PCs[[i]]
  for(j in 1:22){
    curr_cmd = paste("plink2",
                     "--bfile",paste(bfiles,"chr",j,sep=''),
                     "--logistic hide-covar firth-fallback",
                     "--maf 0.01",
                     "--threads 4",
                     paste("--pheno",covar_file),
                     paste("--pheno-name ExerciseGroup"),
                     paste("--covar",covar_file),
                     # ,paste(curr_PCs,collapse=",") # add this below
                     paste("--covar-name sex,Age,",sep=""),
                     "--adjust --out",paste(curr_dir,"chr",j,sep=''))
    run_plink_command(curr_cmd,curr_dir,paste("logistic_chr",j,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000,
                      plink_pkg = "plink/2.0a1")
  }
}
wait_for_job(waittime = 120)

# Merge all results to a single file
out_folders = list.files(gwas_dir)
for(i in 1:out_folders){
  curr_dir = paste(gwas_dir,out_folders[i],"/",sep="")
  curr_files = list.files(curr_dir)
  setwd(curr_dir)
  
  res_files1 = curr_files[grepl("logistic.hybrid$",curr_files)]
  res_file1 = paste(curr_dir,"plink_logistic_assoc.txt",sep="")
  for(j in 1:length(res_files1)){
    if(j==1){
      system(paste("less",res_files1[j],">",res_file1))
    }
    if(j>1){
      system(paste("less",res_files1[j],"| grep -v CHROM >>",res_file1))
    }
  }
  
  res_files2 = curr_files[grepl("adjusted$",curr_files)]
  res_file2 = paste(curr_dir,"plink_logistic_assoc_adj.txt",sep="")
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
for(j in 3:40){
  x = pcax[,j]
  x = abs(x-mean(x))/sd(x)
  outliers_to_rem[x>7] = T
}
tb = table(outliers_to_rem,d$CohortName)
for(j in 1:ncol(tb)){
  print(tb[2,j]/sum(tb[,j]))
}

test_inds = 1:nrow(d)
fisher.test(table(outliers_to_rem[test_inds],
                  d$CohortName[test_inds]=="elite"),alternative = "g")$p.value
fisher.test(table(outliers_to_rem[test_inds],
                  d$CohortName[test_inds]=="cooper"),alternative = "g")$p.value

eu_af_clustering = pcax[,2] > 0.02 
table(eu_af_clustering)
table(eu_af_clustering,d$CohortName)

all_filters = !eu_af_clustering & !rl_clustering & !outliers_to_rem
names(all_filters) = rownames(d)
table(all_filters,d$CohortName)
table(all_filters)
save(eu_af_clustering,rl_clustering,outliers_to_rem,all_filters,file="eu_clustering_elite_data_manual.RData")

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],all_filters[inds],all_filters[inds],
      main = "All filters",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]])

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],eu_af_clustering[inds],eu_af_clustering[inds],
    main = "Non-Asian clustering",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]])

inds = rownames(d)
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],d[inds,]$CohortName,d[inds,]$CohortName,
    main = "By cohort",xlab="PC1",ylab="PC2")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

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

simple_cl = abs(pcax[inds,22]-mean(pcax[inds,22]))/sd(pcax[inds,22]) > 6
res = two_d_plot_visualize_covariate(pcax[inds,22],pcax[inds,22],simple_cl,simple_cl,
    main = "Outliers",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topleft",names(res[[1]]),fill = res[[1]])

par(mfrow=c(2,2))
for(j in seq(1,8,by=2)){
  inds = rownames(d)
  res = two_d_plot_visualize_covariate(pcax[inds,j],pcax[inds,j+1],d[inds,]$CohortName,d[inds,]$CohortName,
      main = "By cohort",xlab=j,ylab=j+1)
  legend(x="topleft",names(res[[1]]),fill = res[[1]],cex=0.5)
  
}
