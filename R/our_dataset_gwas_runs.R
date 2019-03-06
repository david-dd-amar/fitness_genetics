
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/impute2_1000gRef_out/"
chrs = paste("chr",c(1:22),sep="")

# This is the project's raw annotation file
mega_covars_path_raw = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper_with_gp.txt"
# This is a file with imputed sex, computed by our_dataset_prepricessing_flow.R
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"

############################################################################
############################################################################

# These files determine the subject sets and their PCA (i.e., use EU or not)
# All
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/merged_mega_data_autosomal.eigenvec"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_v1/"
system(paste("mkdir",out_path))

# EU
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_v1/"
system(paste("mkdir",out_path))

# EU
pca_results = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/eu_gwas/merged_mega_data_autosomal.eigenvec"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res_eu_v1_wo_age/"
system(paste("mkdir",out_path))


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
mega_covars = mega_covars[rownames(pcax),c(1:5,7:10,14:15,21:24,27)]
mega_covars[mega_covars==""]=NA
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
covs_file = paste(out_path,"pheno.phe",sep="")
write.table(mega_covars,file=covs_file,row.names = F,col.names = T,quote=F,sep=" ")

cooper_pheno_file = paste(out_path,"cooper_pheno.phe",sep="")
write.table(mega_covars[!is.na(cooper_col),c("FID","IID","cooper_col")],file=cooper_pheno_file,row.names = F,col.names = T,quote=F,sep=" ")

elite_pheno_file = paste(out_path,"elite_pheno.phe",sep="")
write.table(mega_covars[!is.na(elite_col),c("FID","IID","elite_col")],file=elite_pheno_file,row.names = F,col.names = T,quote=F,sep=" ")

# # Make sure the covars file fits the fam file
# all(curr_fam[,1]==mega_covars[,1])
# all(curr_fam[,2]==mega_covars[,2])

############################################################################
############################################################################
# Run the GWAS: try different numbers of PCs
covs_file = paste(out_path,"pheno.phe",sep="")
cooper_pheno_file = paste(out_path,"cooper_pheno.phe",sep="")
elite_pheno_file = paste(out_path,"elite_pheno.phe",sep="")

for(num_pcs in 0:7){
  cov_string = "--covar-name sex"
  if(num_pcs > 0){
    cov_string = paste("--covar-name sex,",paste("PC",1:num_pcs,sep="",collapse=","),sep="")
  }
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  system(paste("mkdir",curr_path))
  for(chr in chrs){
    curr_cmd = paste("plink2 --bfile",paste(dataset,chr,sep=''),
                     "--logistic firth-fallback hide-covar",
                     "--pheno",cooper_pheno_file,
                     "--pheno-name cooper_col",
                     "--covar",covs_file,
                     "--maf 0.01",
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
                     "--pheno-name elite_col",
                     "--covar",covs_file,
                     "--maf 0.01",
                     cov_string,
                     "--allow-no-sex --adjust",
                     "--threads",4,
                     "--out",paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chr,sep="")
    )
    run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res_pcs",num_pcs,"_",chr,sep=""),
                      get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000,plink_pkg="plink/2.0a1")
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
for(num_pcs in 0:7){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  res_files = paste(curr_path,"elite_gwas_res_pcs",num_pcs,"_",chrs,".elite_col.glm.logistic.hybrid",sep="")
  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  concatenate_res_files(res_files,res_file)
  res_files = paste(curr_path,"cooper_gwas_res_pcs",num_pcs,"_",chrs,".cooper_col.glm.logistic.hybrid",sep="")
  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
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
for(cc in c("elite","cooper")){
  for(num_pcs in 0:7){
    curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
    all_log_files = list.files(curr_path)
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
for(num_pcs in 0:7){
  curr_path = paste(out_path,"gwas_num_pcs_",num_pcs,"/",sep="")
  
  res_file = paste(curr_path,"elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_elite_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F,comment.char="")
  ps = as.numeric(res[,"P"])
  print(paste("elite 1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
  print(paste("elite 1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
  res = res[,c(1:2,ncol(res))]
  colnames(res) = c("chromosome","position","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  file2pvals[[res_file]] = ps
  
  res_file = paste(curr_path,"cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res_file2 = paste(curr_path,"fuma_cooper_gwas_res_all_pcs",num_pcs,".assoc",sep="")
  res = read.table(res_file,header=T,stringsAsFactors = F,comment.char="")
  ps = as.numeric(res[,"P"])
  print(paste("cooper 1e-8",num_pcs,sum(ps<1e-8,na.rm=T)))
  print(paste("cooper 1e-6",num_pcs,sum(ps<1e-6,na.rm=T)))
  res = res[,c(1:2,ncol(res))]
  colnames(res) = c("chromosome","position","P-value")
  write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")
  file2pvals[[res_file]] = ps
}



