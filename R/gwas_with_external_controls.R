
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
python_script = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/python_sh/recode_indels.py"
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

# Additional files for the PCA analysis
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
sample_metadata_raw = read.delim(sample_metadata,stringsAsFactors = F)
sample_metadata_raw = correct_dups_in_sample_metadata(sample_metadata_raw)

####################################################################################################
####################################################################################################
####################################################################################################
# Create control files
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
wait_for_job()

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

err_path = paste(out_path,"merge_control_beds.err",sep="")
log_path = paste(out_path,"merge_control_beds.log",sep="")
curr_cmd = paste("plink --bfile",all_out_bed_files[1],
                 "--merge-list",allfiles_path,
                 "--make-bed --out",paste(out_path,"merged_control_geno",sep=''))
curr_sh_file = "merged_control_beds.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()
all_out_bed_files = list.files(out_path)
all_out_bed_files = all_out_bed_files[grepl(".bed$",all_out_bed_files)]
system(paste("rm ",out_path,"merge_geno_chr*",sep=""))

# add frequencies
err_path = paste(out_path,"merge_control_beds_frq.err",sep="")
log_path = paste(out_path,"merge_control_beds_frq.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_control_geno",sep=''),
                 "--freq --out",paste(out_path,"merged_control_geno",sep=''))
curr_sh_file = "merged_control_beds_frq.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

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
our_original_snp_ids = as.character(our_bim_data[,2])
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
  # x2 may have redundancies: need to solve this issue
  dups = names(which(table(curr_locs[l1_inds])>1))
  rows_to_remove = c()
  for(dup in dups){
    arr = strsplit(dup,split = ";")[[1]]
    curr_rows = x2[,1]==arr[1] & x2[,4]==arr[2]
    curr_rows = which(curr_rows)[-1]
    rows_to_remove = union(rows_to_remove,curr_rows)
  }
  if(length(rows_to_remove)>0){
    x2 = x2[-rows_to_remove,]  
  }
  x2 = x2[order(as.numeric(x2[,4])),]
  print(all(x1[,4]==x2[,4]))
  our_new_snp_ids[rownames(x1)] = rownames(x2)
}
rownames(our_bim_data) = our_new_snp_ids[rownames(our_bim_data)]
our_bim_data[,2] = rownames(our_bim_data)
table(our_new_snp_ids!=names(our_new_snp_ids)) # shows the number of changed IDs
all(names(our_new_snp_ids)==our_original_snp_ids) # should be TRUE

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
wrong_loc1 = which(intersect_snps_our$V4 != intersect_snps_ukbb$V4)
intersect_snps_ukbb = intersect_snps_ukbb[-wrong_loc1,]
intersect_snps_our = intersect_snps_our[-wrong_loc1,]

############# Correct the Affymetrix indel issue ###############
write.table(t(t(rownames(intersect_snps_our))),
            file=paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)

# 1. Merged geno should use the intersect snps only
err_path = paste(out_path,"merge_control_beds_filt1.err",sep="")
log_path = paste(out_path,"merge_control_beds_filt1.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_control_geno",sep=''),
                 "--extract",paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),
                 "--make-bed --out",paste(out_path,"merged_control_geno",sep=''))
curr_sh_file = "merged_control_beds_filt1.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

# 2. Run our custom script that deals with indels
err_path = paste(out_path,"indels_recode.err",sep="")
log_path = paste(out_path,"indels_recode.log",sep="")
curr_cmd = paste("python",python_script,
  paste(out_path,"merged_control_geno",sep=''),
  paste(out_path,"merged_control_geno",sep='')
)
curr_sh_file = "indels_recode.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

# The result of the script above: the merged_control_geno.bed and its bim
# now code indels using the I/D notation, which fits our Illumina dataset.
# Also, note that the code above will not change a bed that already uses
# the I/D notation.

# For the rest of script below, we now need to reread the new bim.
intersect_snps_ukbb = read.table(paste(out_path,"merged_control_geno.bim",sep=''),
                                 stringsAsFactors = F,header = F)
rownames(intersect_snps_ukbb) = intersect_snps_ukbb[,2]
# sanity checks:
# This may be FALSE due to reordering, but an error is a problem
all(rownames(intersect_snps_ukbb)==rownames(intersect_snps_our))
# These two should be empty:
setdiff(rownames(intersect_snps_ukbb),rownames(intersect_snps_our))
setdiff(rownames(intersect_snps_our),rownames(intersect_snps_ukbb))

# This reformatting is done only for the scripts below,
# it does not mean we assume that the SNP order in both beds is 
# the same
intersect_snps_ukbb = intersect_snps_ukbb[rownames(intersect_snps_our),]

################################################################

# Analysis of the alleles in the two bim files
# 1 means the same alleles
# -1 means the same after reverse complement
# 0 means we cannot explain (hopefully very few of these)
# These can happen, for example, when the same position codes for
# an indel in one platform and a standard SNP in the other platform.
# We ignore these cases. 
# See comment below for an example from the direct genotypes of UKBB.
bim_check = cbind(intersect_snps_our[,5:6],intersect_snps_ukbb[,5:6])
bim_check_res = apply(bim_check,1,check_bims_snp_info)
table(bim_check_res)

# Check those that are zero
# These should be SNPs with the same location in both
# platforms but one is a SNP and the other is an indel.
# These should be ignored.
# An example for such a variant: 
# Affy: Affx-52323799, chr 1, position 57221553
# To recheck this do:
# less /oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/ukb_cal_chr1_v2.bim | grep Affx-52323799
# less /oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/maf_filter.bim | grep 57221553
bim_check[which(bim_check_res==0),]
intersect_snps_our["Affx-52323799",]
intersect_snps_ukbb["Affx-52323799",]

wrong_loc = which(bim_check_res==0)
snps_to_filp = names(bim_check_res)[which(bim_check_res==-1)]

# Another test: did we get indel SNPs to flip?
# All should be false:
table(apply(intersect_snps_our[snps_to_filp,5:6],1,function(x){x[1]=="I"||x[1]=="D"||x[2]=="I"||x[2]=="D"}))

intersect_snps_our = intersect_snps_our[-wrong_loc,]
intersect_snps_ukbb = intersect_snps_ukbb[-wrong_loc,]

# print results into three files: 
#   1. our corrected bim
#     Note that this is the same size of the original bim and in the same 
#     order. The only change is the SNP ids, which can now be compared to
#     the external data.
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

# Below we have a plink-based merging of the data. As an alternative we can analyze merging results
# from other tools such as qctools. These analyses are in external_merge_tools_output_analysis.R script.

# extract the snp intersect and flip snps
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
wait_for_job()

# compare the resulting bims
bimfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_10k_rand_controls_sex_age/our_bed_data_extracted_flipped.bim"
bimfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_10k_rand_controls_sex_age/merged_control_geno.bim"
bim1 = read.table(bimfile1)
bim2 = read.table(bimfile2)
length(intersect(bim1[,2],bim2[,2]))

# merge with our bed file using plink
controls_bed = paste(out_path,"merged_control_geno",sep='')
err_path = paste(out_path,"merge_with_our_bed.err",sep="")
log_path = paste(out_path,"merge_with_our_bed.log",sep="")
curr_cmd = paste("plink --bfile", paste(out_path,"our_bed_data_extracted_flipped",sep=''),
                 "--bmerge",controls_bed,
                 "--extract",paste(out_path,"bims_analysis_intersect_snps.txt",sep=""),
                 "--make-bed --out",paste(out_path,"merged_bed_final_for_gwas",sep=''))
curr_sh_file = "merge_with_our_bed.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

####################################################################################################
####################################################################################################
####################################################################################################
# Run PCA and freq
err_path = paste(out_path,"maf_and_pca.err",sep="")
log_path = paste(out_path,"maf_and_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--pca --freq --out",paste(out_path,"merged_bed_final_for_gwas",sep=''))
curr_sh_file = "maf_and_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

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
table(covars[,"Age"])
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

# 1. Linear of all three groups + sex, age, and 10 PCs
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
err_path = paste(out_path,"gwas_three_groups_linear.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear",sep=''))
curr_sh_file = "gwas_three_groups_linear.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 2. Logistic: Cooper vs. UKBB, + sex, age, and 10 PCs
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
err_path = paste(out_path,"ukbb_vs_cooper.err",sep="")
log_path = paste(out_path,"ukbb_vs_cooper.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--logistic hide-covar firth-fallback",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--1 --adjust",
                 "--out",paste(out_path,"ukbb_vs_cooper",sep=''))
curr_sh_file = "ukbb_vs_cooper.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# 3. Linear of all three groups + sex, and 10 PCs (no age)
pheno_file = paste(out_path,"ukbb_elite_cooper.phe",sep='')
covar_file = paste(out_path,"ukbb_elite_cooper_covars.phe",sep='')
err_path = paste(out_path,"gwas_three_groups_linear_wo_age.err",sep="")
log_path = paste(out_path,"gwas_three_groups_linear_wo_age.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--linear hide-covar",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--adjust",
                 "--out",paste(out_path,"gwas_three_groups_linear_wo_age",sep=''))
curr_sh_file = "gwas_three_groups_linear_wo_age.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()

# 4. Logistic: Elite vs. UKBB, + sex, age, and 10 PCs
sample_inds = as.numeric(covars[,"ExerciseGroup"])!=0
covars_copy = covars
covars_copy[,"ExerciseGroup"] = as.character(as.numeric(covars[,"ExerciseGroup"])+1)
table(covars_copy[sample_inds,"ExerciseGroup"])
pheno_file = paste(out_path,"ukbb_vs_elite.phe",sep='')
covar_file = paste(out_path,"ukbb_vs_elite_covars.phe",sep='')
write.table(file=pheno_file,
            covars_copy[sample_inds,c(1:2,ind)],sep=" ",row.names = F,col.names = T,quote=F)
write.table(file=covar_file,
            covars_copy[sample_inds,-ind],sep=" ",row.names = F,col.names = T,quote=F)
err_path = paste(out_path,"ukbb_vs_elite.err",sep="")
log_path = paste(out_path,"ukbb_vs_elite.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(out_path,"merged_bed_final_for_gwas",sep=''),
                 "--logistic hide-covar firth-fallback",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name ExerciseGroup"),
                 "--allow-no-sex",
                 paste("--covar",covar_file),
                 "--covar-name sex,Batch,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                 "--1 --adjust",
                 "--out",paste(out_path,"ukbb_vs_elite_logistic",sep=''))
curr_sh_file = "ukbb_vs_elite_logistic.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))


# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Repeat the GWAS above after applying PCA and clustering, number of clusters is a parameter
# given by the user. This is determined based on either some prior knowledge or by looking at the
# elbow plot of kmeans.
# PROBLEM (July 2018): PCA only separates UKBB from others
d = read.table(paste(out_path,"ukbb_elite_cooper_covars.phe",sep=''),header=T,stringsAsFactors = F)
d2 = sample_metadata_raw
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
d_cohort = read.table(paste(out_path,"ukbb_elite_cooper.phe",sep=''),header=T,stringsAsFactors = F)

# Association between PCs and UKBB or not
ukbb_inds = d_cohort$ExerciseGroup == -1
pc_pvals = c()
for(k in 1:20){
  x1 = d[ukbb_inds,paste("PC",k,sep="")]
  x2 = d[!ukbb_inds,paste("PC",k,sep="")]
  pc_pvals[k] = wilcox.test(x1,x2)$p.value
}

# Cluster by PCs
# Assumption: use PC 3 and onwards as the first two are
# batch PCs that separate UKBB from non-UKBB
pc_x = as.matrix(d[,c("PC11","PC12")])

# Check number of clusters in PCA plot
wss <- sapply(1:10,
              function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
wss[2:10]/wss[1:9] # by manual inspection we choose 3 clusters here (consistent with 
# the analysis of our samples without UKBB)

num_pca_clusters = 3
pc_x_kmeans = kmeans(pc_x^4,num_pca_clusters)
table(pc_x_kmeans$cluster)
table(pc_x_kmeans$cluster,d_cohort$ExerciseGroup)
kmeans_res = pc_x_kmeans$cluster

# Select subjects from the largest cluster
clustable = table(kmeans_res)
selected_cluster = names(which(clustable == max(clustable)))
selected_subjects = names(which(kmeans_res == selected_cluster))

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# Local analysis
setwd("/Users/David/Desktop/elite/gwas_results/ukbb_pca")
d = read.table("ukbb_elite_cooper_covars.phe",header=T,stringsAsFactors = F)
d2 = read.delim("../../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
table(is_jap)
d_cohort = read.table("ukbb_elite_cooper.phe",header=T,stringsAsFactors = F)
cohorts = d_cohort$ExerciseGroup
cohorts[cohorts=="-1"] = "UKBB"
cohorts[cohorts=="1"] = "ELITE"
cohorts[cohorts=="0"] = "Cooper"
names(cohorts) = d_cohort[,2]
table(cohorts)
all(names(cohorts)==d[,2])

# Cluster by PCs
# Assumption: use PC 3 and onwards as the first two are
# batch PCs that separate UKBB from non-UKBB
pc_x = as.matrix(d[,c("PC11","PC12")])

# Check number of clusters in PCA plot
wss <- sapply(1:10,
              function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
wss[2:10]/wss[1:9] # by manual inspection we choose 3 clusters here (consistent with 
# the analysis of our samples without UKBB)
plot(1:10, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# PCA plots
inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC1",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC2[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC3",ylab="PC2")
legend(x="top",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(d$PC12[inds],
    d$PC11[inds],cohorts[inds],cohorts[inds],
    main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

table(kmeans_res)
res = two_d_plot_visualize_covariate(d$PC12[inds],
    d$PC11[inds],kmeans_res,kmeans_res,
    main = "UKBB, Cooper, ELITE",xlab="PC12",ylab="PC11")
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################

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
                                      paste(out_path,"merged_bed_final_for_gwas.frq",sep=""),maf = 0.001)
  print(dim(res_fuma))
  write.table(res_fuma,
              file= paste(newdir,f,sep=""),
              row.names = F,col.names = T,quote = F,sep=" ")
}

mafs = read.table(paste(out_path,"merged_bed_final_for_gwas.frq",sep=""),
                  stringsAsFactors = F,header=T)
mafs=mafs[,c(1,2,5,6)]
nsample = (mafs[,3]^2 + 2*mafs[,3])*(mafs[,4]/2)
mafs = cbind(mafs,nsample)
write.table(mafs,
            file= paste(newdir,"mafs.txt",sep=""),
            row.names = F,col.names = T,quote = F,sep=" ")




