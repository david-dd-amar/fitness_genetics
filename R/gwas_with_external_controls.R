
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# UKBB: direct genotypes
external_files_path = "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age_with_info.txt"
analysis_name = "ukbb_10k_rand_controls_sex_age"

# UKBB: imputed genotypes
external_files_path = "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age_with_info.txt"
analysis_name = "ukbb_imputed_10k_rand_controls_sex_age"

# Our stuff
our_bed_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/maf_filter"
our_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/three_group_analysis_genepool_controls.phe"
our_phe_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/three_group_analysis_genepool_controls.phe"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/"
chrs = 1:22
all_files = list.files(external_files_path)
if(analysis_name != ""){
  out_path = paste(out_path,analysis_name,"/",sep="")
  system(paste("mkdir",out_path))
}

####################################################################################################
####################################################################################################
####################################################################################################
# Create control files
jobs_before = get_my_jobs()
for (chr in chrs){
  curr_file = all_files[grepl(".bed$",all_files) & grepl(paste("chr",chr,"_",sep=""),all_files)]
  curr_file = gsub(pattern = ".bed$",replacement = "",curr_file)
  err_path = paste(out_path,"merge_geno.err",sep="")
  log_path = paste(out_path,"merge_geno.log",sep="")
  curr_cmd = paste("plink --bfile",paste(external_files_path,curr_file,sep=''),
                   "--keep",external_control_ids,
                   "--make-bed --out",paste(out_path,"merge_geno_chr",chr,sep=""))
  curr_sh_file = "merge_geno.sh"
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}
wait_for_job(jobs_before,5)

####################################################################################################
####################################################################################################
####################################################################################################
# Merge bed/ped control data files
all_out_bed_files = list.files(out_path)
all_out_bed_files = all_out_bed_files[grepl(".bed$",all_out_bed_files)]
all_out_bed_files = gsub(".bed","",all_out_bed_files)
all_out_bed_files = paste(out_path,all_out_bed_files,sep="")
length(all_out_bed_files)
allfiles_path = paste(out_path,"allfiles.txt",sep="")
write.table(t(t(all_out_bed_files[-1])),
            file = allfiles_path,sep="",row.names = F,col.names = F,quote = F)

jobs_before = get_my_jobs()
err_path = paste(out_path,"merge_control_beds.err",sep="")
log_path = paste(out_path,"merge_control_beds.log",sep="")
curr_cmd = paste("plink --bfile",all_out_bed_files[1],
                 "--merge-list",allfiles_path,
                 "--make-bed --out",paste(out_path,"merged_control_geno",sep=''))
curr_sh_file = "merged_control_beds.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=2,mem_size=64000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
list.files(out_path)
all_out_bed_files = list.files(out_path)
all_out_bed_files = all_out_bed_files[grepl(".bed$",all_out_bed_files)]
readLines(log_path)
system(paste("rm ",out_path,"merge_geno_chr*",sep=""))

####################################################################################################
####################################################################################################
####################################################################################################
# analyze the bim files: check SNP intersect, locations,
# and which snps must be flipped before we analyze
# Compare bim files: get SNPs that should be flipped
# TODO: write a custom script that merges the data - will need to recode and solve
# maf issues
our_bim_file = paste(our_bed_path,".bim",sep="")

# Compare bim files: get SNPs that should be flipped
our_bim_data = read.table(our_bim_file,stringsAsFactors = F,header = F)
# correct our bim info
id_is_location = grepl(":",our_bim_data[,2])
# extract true locations from the snp ids
id_is_lc_arr = sapply(our_bim_data[id_is_location,2],function(x)strsplit(x,":|-",perl=T)[[1]][1:2])
our_bim_data[id_is_location,1] = id_is_lc_arr[1,]
our_bim_data[id_is_location,4] = id_is_lc_arr[2,]
our_locs = apply(our_bim_data[,c(1,4)],1,function(x)paste(sort(x),collapse=";"))
rownames(our_bim_data) = our_bim_data[,2]

ukbb_bim_files = list.files(external_files_path)
ukbb_bim_files = ukbb_bim_files[grepl(".bim$",ukbb_bim_files)]

# loop 1: correct our ids to match those in ukbb
our_new_snp_ids = rownames(our_bim_data)
names(our_new_snp_ids) = rownames(our_bim_data)
for(f in ukbb_bim_files){
  curr_bim = read.table(paste(external_files_path,f,sep=""),stringsAsFactors = F)
  rownames(curr_bim) = curr_bim[,2]
  curr_locs = apply(curr_bim[,c(1,4)],1,function(x)paste(sort(x),collapse=";"))
  loc_intersect = intersect(our_locs,curr_locs)
  # analyze the snps with the same location but different ids
  l1_inds = is.element(curr_locs,set=loc_intersect)
  l2_inds = is.element(our_locs,set=loc_intersect)
  x1 = our_bim_data[l2_inds,]
  x1 = x1[order(as.numeric(x1[,4])),]
  x2 = curr_bim[l1_inds,]
  x2 = x2[order(as.numeric(x2[,4])),]
  print(all(x1[,4]==x2[,4]))
  our_new_snp_ids[rownames(x1)] = rownames(x2)
}
rownames(our_bim_data) = our_new_snp_ids[rownames(our_bim_data)]
our_bim_data[,2] = rownames(our_bim_data)
table(our_new_snp_ids!=names(our_new_snp_ids))

# loop2: look at the intersection
intersect_snps_our = c()
intersect_snps_ukbb = c()
for(f in ukbb_bim_files){
  curr_bim = read.table(paste(external_files_path,f,sep=""),stringsAsFactors = F)
  rownames(curr_bim) = curr_bim[,2]
  inds = intersect(rownames(curr_bim),rownames(our_bim_data))
  print(length(inds))
  d1 = our_bim_data[inds,];d2=curr_bim[inds,]
  print(paste(nrow(d1),nrow(d2)))
  intersect_snps_our = rbind(intersect_snps_our,d1)
  intersect_snps_ukbb = rbind(intersect_snps_ukbb,d2)
}
diff_snp_inds = (abs(as.numeric(intersect_snps_our$V4) - as.numeric(intersect_snps_ukbb$V4))>10)
table(diff_snp_inds)
wrong_loc1 = which(intersect_snps_our$V4 != intersect_snps_ukbb$V4)
length(wrong_loc1)

bim_check = cbind(intersect_snps_our[,5:6],intersect_snps_ukbb[,5:6])
bim_check_res = apply(bim_check,1,check_bims_snp_info)
table(bim_check_res)
# bim_check[which(bim_check_res==0),]
# intersect_snps_our["Affx-89021291",]
wrong_loc = union(wrong_loc1,which(bim_check_res==0))
snps_to_filp = names(bim_check_res)[which(bim_check_res==-1)]

intersect_snps_our = intersect_snps_our[-wrong_loc,]
intersect_snps_ukbb = intersect_snps_ukbb[-wrong_loc,]

# print results into three files: 
#   1. our corrected bim
write.table(our_bim_data,file=paste(out_path,"our_bim_data.bim",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)
#   2. Intersect SNPs
write.table(t(t(rownames(intersect_snps_ukbb))),
            file=paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)
#   3. SNPs that should be flipped
write.table(t(t(snps_to_filp)),
            file=paste(out_path,"bims_analysis_intersect_snps_to_flip.txt",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################

# extract the snp intersect and flip snps
jobs_before = get_my_jobs()
err_path = paste(out_path,"extract_snps_and_flip.err",sep="")
log_path = paste(out_path,"extract_snps_and_flip.log",sep="")
curr_cmd = paste("plink --bed",paste(our_bed_path,".bed",sep=""),
                 "--bim", paste(out_path,"our_bim_data.bim",sep=""),
                 "--fam", paste(our_bed_path,".fam",sep=""),
                 "--extract",paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),
                 "--flip", paste(out_path,"bims_analysis_intersect_snps_to_flip.txt",sep=""),
                 "--make-bed --out",paste(out_path,"our_bed_data_extracted_flipped",sep=''))
curr_sh_file = "extract_snps_and_flip.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
readLines(log_path)

# merge with our bed file
controls_bed = paste(out_path,"merged_control_geno",sep='')
jobs_before = get_my_jobs()
err_path = paste(out_path,"merge_with_our_bed.err",sep="")
log_path = paste(out_path,"merge_with_our_bed.log",sep="")
curr_cmd = paste("plink --bfile", paste(out_path,"our_bed_data_extracted_flipped",sep=''),
                 "--bmerge",controls_bed,
                 "--extract",paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),
                 "--make-bed --out",paste(out_path,"merged_bed_final_for_gwas",sep=''))
curr_sh_file = "merge_with_our_bed.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
readLines(log_path)

####################################################################################################
####################################################################################################
####################################################################################################
# Run PCA and freq
jobs_before = get_my_jobs()
err_path = paste(out_path,"maf_and_pca.err",sep="")
log_path = paste(out_path,"maf_and_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--pca --freq --out",paste(out_path,"merged_bed_final_for_gwas",sep=''))
curr_sh_file = "maf_and_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
list.files(out_path)

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
our_phe = as.character(our_covars[,4])
genepool_inds = which(as.numeric(our_phe)==-1)
# all(our_covars[,1]==our_phe[,1])
# table(our_phe[,3])
external_covars = read.table(external_covars_path,stringsAsFactors = F)
external_covars = cbind(as.character(external_covars[,1]),external_covars)
external_samples = as.character(external_covars[,1])

# # read our fam file
# fam_info = read_plink_table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/maf_filter.fam",has_header = F)
# iid_to_fid = fam_info[,1]

# temp sol for ukbb - add batches
batch_data = as.matrix(
  read.table("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/ukb2228.tab.fam",
  stringsAsFactors = F, header=T))
rownames(batch_data) = as.character(batch_data[,1])
external_covars = cbind(batch_data[external_samples,c(1:2,5:6)],external_covars[,4])
external_covars = cbind(external_covars,rep("-1",nrow(external_covars)))
colnames(external_covars) = c("FID","IID","sex","Batch","Age","ExerciseGroup")

# Define the final covariance matrix (with the ExerciseGroup column)
covars = rbind(our_covars[,c(1,2,3,5,6,4)],external_covars)
pca_res = read_plink_table(paste(out_path,"merged_bed_final_for_gwas.eigenvec",sep=''),F)
pca_res = pca_res[,-c(1:2)]
colnames(pca_res) = paste("PC",1:ncol(pca_res),sep='')
table(is.element(covars[,"IID"],set=rownames(pca_res))) # all should be true
covars = cbind(covars,pca_res[covars[,"IID"],])
covars = covars[-genepool_inds,]
covars[,"Batch"] = cov_phe_col_to_plink_numeric_format(covars[,"Batch"])
table(covars[,"Batch"])
for(j in 1:ncol(covars)){
  covars[,j] = gsub(" ","",as.character(covars[,j]))
}

ind = which(colnames(covars)=="ExerciseGroup")
write.table(file=paste(out_path,"ukbb_elite_cooper.phe",sep=''),
            covars[,c(1:2,ind)],sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=paste(out_path,"ukbb_elite_cooper_covars.phe",sep=''),
            covars[,-ind],sep=" ",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# Run GWAS

# 1. Linear of all three groups + sex, age, and 5 PCs
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
#jobs_before = get_my_jobs()
err_path = paste(out_path,"gwas_three_groups_linear.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#list.files(out_path)
#readLines(err_path)

# 2. Logistic: Cooper vs. UKBB, + sex, age, and 5 PCs
sample_inds = as.numeric(covars[,"ExerciseGroup"])!=1
covars_copy = covars
covars_copy[,"ExerciseGroup"] = as.character(as.numeric(covars[,"ExerciseGroup"])+1)
table(covars_copy[sample_inds,"ExerciseGroup"])
pheno_file = paste(out_path,"ukbb_vs_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_vs_cooper_covars.phe",sep='')
write.table(file=pheno_file,
            covars_copy[sample_inds,c(1:2,ind)],sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=covar_file,
            covars_copy[sample_inds,-ind],sep=" ",row.names = F,col.names = T,quote=F)

#jobs_before = get_my_jobs()
err_path = paste(out_path,"gwas_two_groups_logistic.err",sep="")
log_path = paste(out_path,"gwas_two_groups_logistic.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--logistic hide-covar firth-fallback",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5",
                 "--1 --adjust",
                 "--out",paste(out_path,"gwas_two_groups_logistic",sep=''))
curr_sh_file = "gwas_two_groups_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#list.files(out_path)
#readLines(err_path)

# 3. Linear of all three groups + sex, and 5 PCs (no age)
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
#jobs_before = get_my_jobs()
err_path = paste(out_path,"gwas_three_groups_linear_wo_age.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear_wo_age.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear_wo_age",sep=''))
curr_sh_file = "gwas_three_groups_linear_wo_age.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#list.files(out_path)
#readLines(err_path)

# get all result files, transform to fuma and keep in a new directory
res_files = list.files(out_path)
res_files = res_files[grepl("adjusted$",res_files)]
newdir = paste(out_path,"fuma/",sep="")
system(paste("mkdir",newdir))
for (f in res_files){
  res = read.delim(paste(out_path,f,sep=''),stringsAsFactors = F)
  print(table(res$UNADJ < 1e-8))
  res_fuma = from_our_sol_to_fuma_res(paste(out_path,f,sep=""),
                                      paste(out_path,"merged_bed_final_for_gwas.bim",sep=""),
                                      paste(out_path,"merged_bed_final_for_gwas.frq",sep=""),maf = 0.01)
  print(dim(res_fuma))
  write.table(res_fuma,
              file= paste(newdir,f,sep=""),
              row.names = F,col.names = T,quote = F,sep=" ")
}







