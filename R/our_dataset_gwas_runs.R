
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# # Our imputation
# dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/impute2_1000gRef_out/"
# # Michigan server with HRC as reference
# dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/mich_hrc/beds/"
# # Michigan server with 1000G as reference
# dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/mich_1000g/beds/"
# New: august 2019
dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/impute2_1000gRef_out/"

chrs = paste("chr",c(1:22),sep="")

# This is the project's raw annotation file
mega_covars_path_raw = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
# # This is a file with imputed sex, computed by our_dataset_preprocessing_flow.R
# mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"
# This is a file with imputed sex, computed by our_dataset_preprocessing_flow.R
# Updated: July 2019
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/integrated_sample_metadata_and_covariates.phe"

# These files determine the subject sets and their PCA (i.e., use EU or not)
# Also, the output path

# # All
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_v1/"
# system(paste("mkdir",out_path))
# 
# # EU
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_v1/"
# system(paste("mkdir",out_path))
# 
# # EU WO age
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_v1_wo_age/"
# system(paste("mkdir",out_path))
# 
# # EU + Michigan impute + HRC
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_mich_hrc/"
# system(paste("mkdir",out_path))
# 
# # EU + Michigan impute + HRC
# pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_mich_1000g/"
# system(paste("mkdir",out_path))

# New run: August 2019
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/eu_gwas/merged_mega_data_autosomal.eigenvec"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/eu_gwas/"
system(paste("mkdir",out_path))
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/eu_gwas_maf0.05/"
system(paste("mkdir",out_path))

# Define the pheno and cov files
covs_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/eu_gwas/all_covs_and_pheno.phe"

############################################################################
############################################################################
############################################################################
############################################################################
cooper_pheno_file = covs_file
elite_pheno_file = covs_file

# Run the GWAS: try different numbers of PCs
for(num_pcs in 0:7){
  cov_string = "--covar-name sex,age"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  system(paste("mkdir",curr_path))
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",cooper_pheno_file,
                     "--pheno-name cooper_vs_gp",
                     "--covar",covs_file,
                     "--maf 0.05",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
  }
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",elite_pheno_file,
                     "--pheno-name elite_vs_gp",
                     "--covar",covs_file,
                     "--maf 0.05",
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
                     "--maf 0.05",
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
                     "--maf 0.05",
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
                     "--maf 0.05",
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
############################################################################
# Analyze the results

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

gwas_job_names = c("cooper_gwas_res_pcs",
              "elite_gwas_res_pcs",
              "elite_vo2_ml_kg_min_PCs",
              "elite_vo2max_l_PCs",
              "cooper_treadmill_time_PCs")

for(num_pcs in 0:7){
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
  for(num_pcs in 0:7){
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
for(num_pcs in 7:0){
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
    write.table(fuma_res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  }
}



