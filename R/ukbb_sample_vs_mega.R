
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

maf_filter_line = "--maf 0.05"
geno_filter_line = "--geno 0.01"
ld_filter_line = paste("--indep-pairwise 500 10",0.3)
tested_pcs_sets = c(5,10,20)

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

# Get the palindromic snps, this may be a source of bias
is_snp_paly<-function(x){
  return(x=="AT" || x=="TA" || x=="CG" || x=="GC")
}

pali_snps_ids = c()
for(chr in chrs){
  # Get the palindromic SNPs
  currbim = read.table(paste(bfiles,"chr",chr,".bim",sep=''),stringsAsFactors = F)
  alleles2 = paste(currbim[,5],currbim[,6],sep="")
  pali_snps2 = sapply(alleles2,is_snp_paly)
  pali_snps2 = currbim[pali_snps2,2]
  pali_snps2 = as.character(pali_snps2)
  pali_snps_ids = c(pali_snps_ids,pali_snps2)
  print(length(pali_snps_ids))
}
pali_snps_file = paste(bfiles,"all_pali_snps.txt",sep="")
write.table(t(t(pali_snps_ids)),file=pali_snps_file,row.names = F,col.names = F,quote = F)

# to include pali snps, change to ""
exclude_pali_plink_line = paste("--exclude",pali_snps_file)

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
                   maf_filter_line,
                   geno_filter_line,
                   ld_filter_line,
                   exclude_pali_plink_line,
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

length(readLines(paste(ld_prune_res,"merged_dataset.bim",sep='')))

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
eu_samples = get(load(paste(out_path,"../ld_prune/manual_clustering.RData",sep="")))
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

# Run the GP vs. UKBB GWAS: try different numbers of PCs
chrs= paste("chr",1:22,sep="")
for(num_pcs in 0:4){
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
                     maf_filter_line,
                     geno_filter_line,
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
for(num_pcs in 0:4){
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

# March 2019: taked the 3 PCs results
for(nn in c("3")){
  
  curr_dir = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/",sep="")
  system(paste("mkdir",curr_dir))
  
  surr_snps_to_exclude = gp_vs_ukbb_res[[nn]]
  # if we wish to exclude pali snps, add them to the list
  # of excluded snps
  if(exclude_pali_plink_line!=""){
    pali_snps_ids = read.table(pali_snps_file,stringsAsFactors = F)
    surr_snps_to_exclude = c(surr_snps_to_exclude,pali_snps_ids[,1])
  }
  print(paste("total excluded snps:",length(surr_snps_to_exclude)))
  curr_snps_to_exclude_file = paste(curr_dir,"surr_snps_to_exclude.txt",sep="")
  write.table(t(t(surr_snps_to_exclude)),file=curr_snps_to_exclude_file,row.names = F,col.names = F,
              quote = F)
  
  # LD prune, PCA on all samples
  ld_prune_res = paste(curr_dir,"ld_prune/",sep="")
  system(paste("mkdir",ld_prune_res))
  
  # Run LD prunning on each bed file
  for(chr in chrs){
    curr_name = paste("ld_",chr,sep='')
    curr_cmd = paste("plink --bfile",paste(bfiles,"chr",chr,sep=''),
                     "--threads 8",
                     maf_filter_line,
                     geno_filter_line,
                     ld_filter_line,
                     "--keep", curr_sample_file,
                     "--exclude",curr_snps_to_exclude_file,
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
                     "--exclude",curr_snps_to_exclude_file,
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

# April 2019: do not exclude the pali snps as in the previous
# analyses
for(nn in c("3")){
  curr_dir = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/",sep="")
  curr_snps_to_exclude_file = paste(curr_dir,"surr_snps_to_exclude.txt",sep="")
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
  for(num_pcs in tested_pcs_sets){
    cov_string = paste("--covar-name sex,Age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")

    # elite
    for(chr in chrs){
        curr_cmd = paste("plink2 --bfile",paste(bfiles,chr,sep=''),
                       "--logistic firth-fallback hide-covar",
                       "--pheno",curr_cov_file,
                       "--pheno-name elite_vs_ukbb",
                       "--covar",curr_cov_file,
                       maf_filter_line,
                       geno_filter_line,
                       cov_string,
                       "--allow-no-sex --adjust",
                       "--keep", curr_sample_file,
                       # "--exclude",curr_snps_to_exclude_file,
                       "--threads",4,
                       "--out",paste(curr_dir,"elite_vs_ukbb_pcs",num_pcs,"_",chr,sep="")
        )
        run_plink_command(curr_cmd,curr_dir,paste("elite_vs_ukbb",num_pcs,"_",chr,sep=""),
                          get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
    }
    
    # cooper
    for(chr in chrs){
      curr_cmd = paste("plink2 --bfile",paste(bfiles,chr,sep=''),
                       "--logistic firth-fallback hide-covar",
                       "--pheno",curr_cov_file,
                       "--pheno-name cooper_vs_ukbb",
                       "--covar",curr_cov_file,
                       maf_filter_line,
                       geno_filter_line,
                       cov_string,
                       "--allow-no-sex --adjust",
                       "--keep", curr_sample_file,
                       # "--exclude",curr_snps_to_exclude_file,
                       "--threads",4,
                       "--out",paste(curr_dir,"cooper_vs_ukbb_pcs",num_pcs,"_",chr,sep="")
      )
      run_plink_command(curr_cmd,curr_dir,paste("cooper_vs_ukbb",num_pcs,"_",chr,sep=""),
                        get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
    }
    
    # GP
    for(chr in chrs){
      curr_cmd = paste("plink2 --bfile",paste(bfiles,chr,sep=''),
                       "--logistic firth-fallback hide-covar",
                       "--pheno",curr_cov_file,
                       "--pheno-name gp_vs_ukbb",
                       "--covar",curr_cov_file,
                       maf_filter_line,
                       geno_filter_line,
                       cov_string,
                       "--allow-no-sex --adjust",
                       "--keep", curr_sample_file,
                       # "--exclude",curr_snps_to_exclude_file,
                       "--threads",4,
                       "--out",paste(curr_dir,"gp_vs_ukbb_pcs",num_pcs,"_",chr,sep="")
      )
      run_plink_command(curr_cmd,curr_dir,paste("gp_vs_ukbb",num_pcs,"_",chr,sep=""),
                        get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
    }
  }
}

####################################################################################################
####################################################################################################
####################################################################################################
# Create FUMA files

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

# March 2019: we take the GP vs. UKBB adjusted for 3 PCs
nn = "3"

for(num_pcs in tested_pcs_sets){
  curr_path = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
  all_files = list.files(curr_path)
  all_files = all_files[grepl(paste("pcs",num_pcs,sep=""),all_files) &
                          grepl("elite_vs_ukbb.glm.logistic.hybrid$",all_files)]
  res_files = paste(curr_path,all_files,sep="")
  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
  all_files = list.files(curr_path)
  all_files = all_files[grepl(paste("pcs",num_pcs,sep=""),all_files) &
                          grepl("cooper_vs_ukbb.glm.logistic.hybrid$",all_files)]
  res_files = paste(curr_path,all_files,sep="")
  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
  all_files = list.files(curr_path)
  all_files = all_files[grepl(paste("pcs",num_pcs,sep=""),all_files) &
                          grepl("gp_vs_ukbb.glm.logistic.hybrid$",all_files)]
  res_files = paste(curr_path,all_files,sep="")
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

all_lambdas = list()
all_lambdas[["elite"]] = list()
all_lambdas[["cooper"]] = list()
for(cc in c("elite","cooper","gp")){
  for(num_pcs in tested_pcs_sets){
    curr_path = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
    all_log_files = list.files(curr_path)
    all_log_files = all_log_files[grepl(paste("pcs",num_pcs,sep=""),all_log_files)]
    all_log_files = all_log_files[grepl("log$",all_log_files)]
    all_log_files = all_log_files[grepl("pcs",all_log_files)]
    curr_files = all_log_files[grepl(num_pcs,all_log_files) & grepl(cc,all_log_files)]
    curr_files = paste(curr_path,curr_files,sep="")
    curr_lambdas = sapply(curr_files, get_lambda_from_log_file)
    print(paste(cc,num_pcs,mean(curr_lambdas,na.rm=T)))
    all_lambdas[[cc]][[as.character(num_pcs)]]=curr_lambdas
  }
}
sapply(all_lambdas,sapply,median,na.rm=T)
save(all_lambdas,file=paste(out_path,"all_lambdas.RData",sep=""))

# reformat and write to fuma files
file2pvals = c()
mhc_region = c(25000000,35000000)
for(num_pcs in tested_pcs_sets){
  curr_path = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F,comment.char="")
  # exclude mhc
  is_mhc = res[,1]==6 & res[,2]>mhc_region[1] & res[,2]<mhc_region[2]
  res = res[!is_mhc,]
  ps = as.numeric(res[,"P"])
  print(paste("elite 1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
  print(paste("elite 1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
  res = res[,c(1:2,ncol(res))]
  colnames(res) = c("chromosome","position","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  file2pvals[[res_file]] = res
  
  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res=NULL
  try({res = read.table(res_file,header=T,stringsAsFactors = F,comment.char="")})
  if(is.null(res)){res = read.delim(res_file,header=T,stringsAsFactors = F,comment.char="")}
  # exclude mhc
  is_mhc = res[,1]==6 & res[,2]>mhc_region[1] & res[,2]<mhc_region[2]
  res = res[!is_mhc,]
  ps = as.numeric(res[,"P"])
  print(paste("cooper 1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
  print(paste("cooper 1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
  res = res[,c(1:2,ncol(res))]
  colnames(res) = c("chromosome","position","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  file2pvals[[res_file]] = res
  
  res_file = paste(curr_path,"gp_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_gp_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res=NULL
  try({res = read.table(res_file,header=T,stringsAsFactors = F,comment.char="")})
  if(is.null(res)){res = read.delim(res_file,header=T,stringsAsFactors = F,comment.char="")}
  # exclude mhc
  is_mhc = res[,1]==6 & res[,2]>mhc_region[1] & res[,2]<mhc_region[2]
  res = res[!is_mhc,]
  ps = as.numeric(res[,"P"])
  print(paste("GP 1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
  print(paste("GP 1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
  res = res[,c(1:2,ncol(res))]
  colnames(res) = c("chromosome","position","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  file2pvals[[res_file]] = res
}

# Compare two results
names(file2pvals)[8:9]
x1 = file2pvals[[8]] 
x2 = file2pvals[[9]]
all(x1[,2]==x2[,2])
p1 = as.numeric(x1[,3])
p2 = as.numeric(x2[,3])
inds = !is.na(p1) & !is.na(p2)
cor(p1[inds],p2[inds])
table(p1[inds]<1e-8 & p2[inds] > 0.01)
table(p2[inds]<1e-8 & p1[inds] > 0.01)

cor.test(-log(p1[inds]),-log(p2[inds]))
cor(-log(p1[inds]),-log(p2[inds]))

curr_path = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
res_file2 = paste(curr_path,"fuma_cooper_gwas_res_20pcs_cleaned_1e6_0.01.assoc",sep="")
inds2 = p1 < 1e-6 & p2 > 0.01
inds2[is.na(inds2)]=F
sum(inds2,na.rm=T)
fuma_res = x1[inds2,]
write.table(fuma_res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
# write.table(fuma_res,col.names = T,row.names = F,quote = F,sep=" ")

# Compare two results
x1 = file2pvals[[5]] 
x2 = file2pvals[[6]]
all(x1[,2]==x2[,2])
p1 = as.numeric(x1[,3])
p2 = as.numeric(x2[,3])
inds = !is.na(p1) & !is.na(p2)
curr_path = paste(bfiles,"eu_gp_vs_ukbb_analysis_",nn,"/gwas/",sep="")
res_file2 = paste(curr_path,"fuma_cooper_gwas_res_10pcs_cleaned_1e6_0.01.assoc",sep="")
inds2 = p1 < 1e-6 & p2 > 0.01
inds2[is.na(inds2)]=F
fuma_res = x1[inds2,]
write.table(fuma_res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")

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
# pcax = read_pca_res("../mega_ukbb_eu/run_pca.eigenvec") # old results for reference

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
            xlab="PC1",ylab="PC2",lwd=3,cex=curr_cex)
legend(x="topleft",names(res[[1]]),fill = res[[1]])
par(mfrow=c(2,2))
res = two_d_plot_visualize_covariate(pcax[inds,3],pcax[inds,4],curr_cohorts[inds],curr_cohorts[inds],
            xlab="PC3",ylab="PC4",lwd=3,cex=curr_cex)
legend(x="topleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,5],pcax[inds,6],curr_cohorts[inds],curr_cohorts[inds],
        xlab="PC5",ylab="PC6",lwd=3,cex=curr_cex)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,7],pcax[inds,8],curr_cohorts[inds],curr_cohorts[inds],
        xlab="PC7",ylab="PC8",lwd=3,cex=curr_cex)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,9],pcax[inds,10],curr_cohorts[inds],curr_cohorts[inds],
        xlab="PC9",ylab="PC10",lwd=3,cex=curr_cex)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,11],pcax[inds,12],curr_cohorts[inds],curr_cohorts[inds],
        xlab="PC9",ylab="PC10",lwd=3,cex=curr_cex)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])
res = two_d_plot_visualize_covariate(pcax[inds,13],pcax[inds,14],curr_cohorts[inds],curr_cohorts[inds],
        xlab="PC9",ylab="PC10",lwd=3,cex=curr_cex)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# plot a single pc
pc=1
par(mfrow=c(2,1))
hist(pcax[inds[curr_cohorts[inds]=="UKBB"],pc],breaks=200,main="PC1 UKBB samples")
hist(pcax[inds[curr_cohorts[inds]!="UKBB"],pc],breaks=200, main = "PC1 non-UKBB samples")

# Do the pcs perfectly predict the group?
y = as.factor(curr_cohorts=="ukbb")
x = pcax[,1:20]
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
manual_clustering = pcax[,1]< 0.01 & pcax[,2] > -0.02 & pcax[,2]<0.02
res = two_d_plot_visualize_covariate(
  pcax[inds,1],pcax[inds,2],manual_clustering[inds],manual_clustering[inds],
  xlab="PC1",ylab="PC2",lwd=3,cex=curr_cex)
table(curr_cohorts[manual_clustering])
save(manual_clustering,file="manual_clustering.RData")

# manual_clustering2 = pcax[,1]< 0.01 & pcax[,2] > -0.02 & pcax[,2]<0.02
# table(manual_clustering,manual_clustering2)

# compare to our mega-only clustering
xx1 = get(load("/Users/David/Desktop/elite/november2018_analysis/mega_with_gp/manual_clustering.RData"))
xx2 = manual_clustering[names(xx1)]
table(xx1,xx2)





