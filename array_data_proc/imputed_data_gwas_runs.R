
chrs = paste("chr",c(1:22),sep="")
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# # Our imputation
# dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/impute2_1000gRef_out/"
# New: august 2019
dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/impute2_1000gRef_out/"

# This is the project's raw annotation file
mega_covars_path_raw = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
# # This is a file with imputed sex, computed by our_dataset_preprocessing_flow.R
# mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"
# This is a file with imputed sex, computed by our_dataset_preprocessing_flow.R
# Updated: July 2019
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/integrated_sample_metadata_and_covariates.phe"

# These files determine the subject sets and their PCA (i.e., use EU or not)
# Also, the output path.

# # EU
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_v1/"
# system(paste("mkdir",out_path))
# 
# New run: August 2019
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/eu_gwas/merged_mega_data_autosomal.eigenvec"
# 5% maf
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/eu_gwas/"
system(paste("mkdir",out_path))
# 1% maf
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/eu_gwas_maf0.01/"
system(paste("mkdir",out_path))

# Define the pheno and cov files
covs_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/eu_gwas/all_covs_and_pheno.phe"

# # Add Cooper vs. ELITE
# covs = read.table(covs_file,header = T,sep=" ",stringsAsFactors = F)
# elite_vs_cooper = covs$elite_vs_gp
# cooper_inds = is.na(elite_vs_cooper)
# elite_vs_cooper[elite_vs_cooper==1] = NA
# elite_vs_cooper[cooper_inds] = 1
# table(elite_vs_cooper)
# covs$elite_vs_cooper = elite_vs_cooper
# write.table(covs,file=covs_file,sep=" ",quote = F,row.names = F,col.names = T)

tested_pcs = c(5,10,20)
maf_line = "--maf 0.01"

############################################################################
############################################################################
############################################################################
cooper_pheno_file = covs_file
elite_pheno_file = covs_file

# Run the GWAS: try different numbers of PCs
for(num_pcs in tested_pcs){
  cov_string = "--covar-name sex,age"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  system(paste("mkdir",curr_path))
  
  # Elite vs. Cooper
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",elite_pheno_file,
                     "--pheno-name elite_vs_cooper",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_vs_cooper_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_vs_cooper_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  
  # Cooper vs. GP
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",cooper_pheno_file,
                     "--pheno-name cooper_vs_gp",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  
  # Elite vs. GP
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",elite_pheno_file,
                     "--pheno-name elite_vs_gp",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  
  # VO2
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--linear hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name VO2max..ml.kg.min.",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_vo2_ml_kg_min_PCs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_vo2_ml_kg_min_PCs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--linear hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name VO2max..l.",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_vo2max_l_PCs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_vo2max_l_PCs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  # Run test
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--linear hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name Treadmill.time",
                     "--covar",covs_file,
                     maf_line,
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"cooper_treadmill_time_PCs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("cooper_treadmill_time_PCs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
}

############################################################################
############################################################################
############################################################################
# Analyze the results

# Merge the results into single files
concatenate_res_files<-function(res_files,res_file){
  for(j in 1:length(res_files)){
    if(j==1){
      system(paste("less",res_files[j],">",res_file))
    }
    if(j>1){
      system(paste("less",res_files[j],"| grep -v CH >>",res_file))
    }
  }
}

gwas_job_names = c(
              "elite_vs_cooper_gwas_res_pcs",
              "cooper_gwas_res_pcs",
              "elite_gwas_res_pcs",
              "elite_vo2_ml_kg_min_PCs",
              "elite_vo2max_l_PCs",
              "cooper_treadmill_time_PCs"
)

for(num_pcs in tested_pcs){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  all_path_files = list.files(curr_path)
  for(jobname in gwas_job_names){
    res_files = all_path_files[grepl(jobname,all_path_files)]
    res_files = res_files[grepl("glm",res_files)]
    res_files = res_files[!grepl("\\.id$",res_files)]
    res_files = res_files[!grepl("adjusted$",res_files)]
    res_files = paste(curr_path,res_files,sep="")
    print(res_files)
    res_file = paste(curr_path,jobname,num_pcs,"_all.assoc",sep="")
    concatenate_res_files(res_files,res_file)
  }
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
for(cc in gwas_job_names){
  all_lambdas[[cc]] = list()
  for(num_pcs in tested_pcs){
    curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
    all_log_files = list.files(curr_path)
    all_log_files = all_log_files[grepl("log$",all_log_files)]
    all_log_files = all_log_files[grepl("pcs",all_log_files,ignore.case = T)]
    curr_files = all_log_files[grepl(num_pcs,all_log_files) & grepl(cc,all_log_files)]
    curr_files = paste(curr_path,curr_files,sep="")
    # print(curr_files);break
    curr_lambdas = sapply(curr_files, get_lambda_from_log_file)
    print(paste(cc,num_pcs,mean(curr_lambdas,na.rm=T)))
    all_lambdas[[cc]][[as.character(num_pcs)]]=curr_lambdas
  }
}
sapply(all_lambdas,sapply,median,na.rm=T)
save(all_lambdas,file=paste(out_path,"all_lambdas.RData",sep=""))

# reformat and write to fuma files
library(data.table,lib.loc = "~/R/packages/")
for(num_pcs in tested_pcs){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  for(jobname in gwas_job_names){
    res_file = paste(curr_path,jobname,num_pcs,"_all.assoc",sep="")
    res_file2 = paste(curr_path,"fuma_",jobname,num_pcs,"_all.assoc",sep="")
    res = fread(res_file,header=T,stringsAsFactors = F,data.table = F)
    ps = as.numeric(res[,"P"])
    print(paste(jobname,"1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
    print(paste(jobname,"1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
    colnames(res)[1:2] = c("chromosome","position")
    colnames(res)[colnames(res)=="P"] = "P-value"
    fuma_res = res
    fwrite(fuma_res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  }
}

# ELITE: meta-analysis using plink
meta_analysis_jobs = c(
  "elite_vs_cooper_gwas_res_pcs",
  "elite_gwas_res_pcs"
)

# Use plink1
library(data.table,lib.loc = "~/R/packages/")
for(num_pcs in tested_pcs){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  # get the results files
  files_str = ""
  for(jobname in meta_analysis_jobs){
    res_file = paste(curr_path,jobname,num_pcs,"_all.assoc",sep="")
    # print an alternative for plink1
    res = fread(res_file,header=T,stringsAsFactors = F,data.table = F)
    colnames(res)[3] = "SNP"
    colnames(res)[1] = "CHR"
    colnames(res)[2] = "BP"
    colnames(res)[4] = "A2"
    res_file2 = paste(curr_path,jobname,num_pcs,"_all_plink1.assoc",sep="")
    files_str = paste(files_str,res_file2,sep=" ")
    fwrite(res,file=res_file2,sep="\t",quote = F,row.names = F,col.names = T)
  }
  curr_cmd = paste("plink --meta-analysis",files_str,
                   "--threads",8,
                   "--out",paste(curr_path,"elite_metaanalysis_PCs",num_pcs,sep="")
  )
  run_plink_command(curr_cmd,curr_path,paste("elite_metaanalysis_PCs",num_pcs,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
}

# Read the results of the meta-analysis
for(num_pcs in tested_pcs){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  res_file = paste(curr_path,"elite_metaanalysis_PCs",num_pcs,".meta",sep="")
  res = fread(res_file,header=T,stringsAsFactors = F,data.table = F)
  print(table(res$P < 5e-08))
  res_file2 = paste(curr_path,"fuma_elite_metaanalysis_PCs",num_pcs,".meta",sep="")
  colnames(res)[1:2] = c("chromosome","position")
  colnames(res)[colnames(res)=="P"] = "P-value"
  fwrite(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
}

# # ELITE: meta-analysis in R
# library(metafor,lib.loc = "~/R/packages/")
# library(data.table,lib.loc = "~/R/packages/")
# library(parallel)
# meta_analysis_jobs = c(
#   "elite_vs_cooper_gwas_res_pcs",
#   "elite_gwas_res_pcs")
# run_simple_meta <- function(x){
#   n = length(x)
#   yi = x[1:(n/2)]
#   vi = x[((n/2)+1):n]
#   m = NULL
#   try({m = rma.uni(yi,vi)})
#   if(is.null(m)){
#     v = c(0,1e-05,1)
#   }
#   else{
#     v = c(m$beta[1],m$se[1],m$pval[1])
#   }
#   names(v) = c("log_OR","se","pval")
#   return(v)
# }
for(num_pcs in tested_pcs){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  Yi = c()
  Vi = c()
  P = c()
  Z = c()
  snp_names = c()
  for(jobname in meta_analysis_jobs){
    res_file = paste(curr_path,jobname,num_pcs,"_all.assoc",sep="")
    res = fread(res_file,header=T,stringsAsFactors = F,data.table = F)
    stats = log(as.numeric(res$OR))
    stats_var = (as.numeric(res$SE))^2
    names(stats) = res[,3]
    names(stats_var) = names(stats)
    Yi = cbind(Yi,stats)
    Vi = cbind(Vi,stats_var)
    P = cbind(P,as.numeric(res$P))
    Z = cbind(Z,as.numeric(res$Z_STAT))
    snp_names = cbind(snp_names,res$ID)
  }
  # QA
  # table(apply(snp_names, 1, function(x)all(x==x[1])))
  meta_stats = cbind(Yi,Vi)
  # ignore these rows
  to_rem = apply(is.na(meta_stats),1,any) | (rowSums(Vi) == 0)
  # run the analysis
  cl = makeCluster(detectCores(),type="FORK")
  system.time({meta_results = t(parLapply(cl,meta_stats[!to_rem,],run_simple_meta))})
  save(meta_results,file=paste(curr_path,"elite_meta_results.RData",sep=""))
}








