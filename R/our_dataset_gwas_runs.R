
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/impute2_1000gRef_out/"
chrs = paste("chr",c(1:22),sep="")
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_v1/"
system(paste("mkdir",out_path))

# This is the project's raw annotation file
mega_covars_path_raw = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper_with_gp.txt"
# This is a file with imputed sex, computed by our_dataset_prepricessing_flow.R
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"


############################################################################
############################################################################

# These files determine the subject sets and their PCA (i.e., use EU or not)
# All
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/merged_mega_data_autosomal.eigenvec"
# EU
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"


############################################################################
############################################################################
# Create the phe files
mega_anno = read.delim(mega_covars_path_raw,stringsAsFactors = F) # main annotation file
mega_anno = mega_anno[!is.na(mega_anno$Sample_ID),]
tmp = paste(mega_anno[,1],mega_anno[,2],sep="_")
mega_covars = read.table(mega_covars_path,stringsAsFactors = F,header=T) # for imputed sex from geno data
curr_fam = read.table(paste(dataset,"chr1.fam",sep=''))
rownames(curr_fam) = curr_fam$IID
rownames(mega_covars) = mega_covars$IID

# remove duplications in the raw annotation file
to_rem = rep(F,nrow(mega_anno))
for(id in names(which(table(tmp)>1))){
  curr_inds = which(tmp==id)
  curr_m = mega_anno[curr_inds,]
  curr_ind = curr_inds[is.element(curr_m$Sample_ID,set=mega_covars$Sample_ID)]
  curr_inds = setdiff(curr_inds,curr_ind)
  to_rem[curr_inds]=T
}
mega_anno = mega_anno[!to_rem,]
rownames(mega_anno) = tmp[!to_rem]

# Update the mega_covars based on the raw data
gp_samples = rownames(mega_covars)[mega_covars$Cohort == "genepool"]
gp_samples = gp_samples[!is.na(gp_samples)]
mega_covars[gp_samples,"age"] = as.numeric(mega_anno[gp_samples,"Age..at.test."])
mega_covars = mega_covars[,!grepl("^PC",colnames(mega_covars))]
colnames(mega_covars)[grepl("Ethni",colnames(mega_covars))] = "Ethnicity"
mega_covars[gp_samples,"Ethnicity"] = mega_anno[gp_samples,"Ethnicity"]
pcax = read_pca_res(pca_results)
mega_covars = mega_covars[rownames(pcax),]
mega_covars = cbind(mega_covars,pcax)


# Define the subcohorts and the pheno files
curr_cohorts = mega_covars$Cohort
cooper_col = rep(NA,nrow(mega_covars))
names(cooper_col) = rownames(mega_covars)
cooper_col[curr_cohorts=="1"]  = 1
cooper_col[curr_cohorts !="1" & curr_cohorts!= "2"]  = 2
elite_col = rep(NA,nrow(mega_covars))
names(elite_col) = rownames(mega_covars)
elite_col[curr_cohorts=="2"]  = 1
elite_col[curr_cohorts !="1" & curr_cohorts!= "2"]  = 2
mega_covars = cbind(mega_covars,elite_col,cooper_col)

# write the pheno file
covs_file = paste(out_path,"pheno.phe")
write.table(mega_covars,file=covs_file,row.names = F,col.names = T,quote=F,sep=" ")

############################################################################
############################################################################
# Run the GWAS: try different numbers of PCs
for(num_pcs in 0:5){
  cov_string = "--covar-name sex,age"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,age,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  system(paste("mkdir",curr_path))
  for(chr in chrs){
    curr_cmd = paste("plink --bfile",paste(dataset,chr,sep=''),
                     "--logistic hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name cooper_col",
                     "--covar",covs_file,
                     "--maf 0.05",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
  }
  for(chr in chrs){
    curr_cmd = paste("plink --bfile",paste(dataset,chr,sep=''),
                     "--logistic hide-covar",
                     "--pheno",covs_file,
                     "--pheno-name elite_col",
                     "--covar",covs_file,
                     "--maf 0.05",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res_pcs",num_pcs,"_",chr,sep=""),
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

