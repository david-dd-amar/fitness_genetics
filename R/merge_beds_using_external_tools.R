# In this script we take two bed files and merge them as follows:
# Input is two bed paths that are assumed to be comparable (no strand and MAF issues)
# 1. We use plink to transform the files to bgen
# 2. We merge the bgens using qctool
# 3. We convert the bgen to bed - the output
# 4. (Optional): given a case-control phenotype we check for issues using PLINK's flip scan
# As for the input, note that the analysis is not symmetric. Whenever possible, we convert the ids
# in file2 to match those in file1. For example if we use ukbb for file1 and Illumina ids for 
# file2 then non-standard ids such as (examXXX) will be mapped to ids from ukbb.

# # Original analysis
# bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-updated"
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/final_dataset_for_analysis-updated"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/"
# 
# # Elite only
# bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-updated"
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/our_prepro/final_dataset_for_analysis-updated"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/with_ukbb/"

# September 2018: new MEGA analysis, HRC as the panel
bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-updated"
bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data_autosomal-hrc_updated"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/"

try(system(paste("mkdir",out_path)))

check_bim_info = T
qctool_path = "/home/users/davidama/apps/qctool_v2/build/release/qctool_v2.0.1"

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

####################################################################################################
####################################################################################################
####################################################################################################
# Compare the bim files
# (1) Check SNP intersect, locations, and which snps must be flipped before we analyze
process_bim_data<-function(bfile1){
  bim_data1 = read.table(paste(bfile1,".bim",sep=""),stringsAsFactors = F,header = F)
  bim_data1_snp_ids = as.character(bim_data1[,2])
  # correct our bim info if needed
  id_is_location = grepl(":",bim_data1[,2])
  print(paste("num snp ids that are location:",sum(id_is_location)))
  # extract true locations from the snp ids
  id_is_lc_arr = sapply(bim_data1[id_is_location,2],function(x)strsplit(x,":|-",perl=T)[[1]][1:2])
  bim_data1[id_is_location,1] = id_is_lc_arr[1,]
  bim_data1[id_is_location,4] = id_is_lc_arr[2,]
  rownames(bim_data1) = bim_data1[,2]
  return(list(bim_data1,id_is_location))
}
bim_data1 = process_bim_data(bfile1)
bim_data2 = process_bim_data(bfile2)
shared_snps = intersect(rownames(bim_data2[[1]]),rownames(bim_data1[[1]]))
intersected_locations = list()
for(chr in unique(bim_data2[[1]][,1])){
  rows1 = bim_data1[[1]][,1]==chr
  rows2 = bim_data2[[1]][,1]==chr
  x1 = bim_data1[[1]][rows1,]
  x2 = bim_data2[[1]][rows2,]
  curr_loc_intersect = intersect(x1[,4],x2[,4])
  curr_snp_id_intersect = intersect(x1[,2],x2[,2])
  print(paste(length(curr_loc_intersect),length(curr_snp_id_intersect)))
  y2_1 = is.element(x2[,4],set=curr_loc_intersect) 
  y2_2 = is.element(x2[,2],set=curr_snp_id_intersect)
  positions_in_intersect_no_id = x2[y2_1&!y2_2,4]
  rownames(x1) = x1[,4]
  rownames(x2) = x2[,4]
  x1 = x1[positions_in_intersect_no_id,]
  x2 = x2[positions_in_intersect_no_id,]
  intersected_locations[[chr]] = list(curr_snp_id_intersect,positions_in_intersect_no_id,x1,x2)
}
save(intersected_locations,file=paste(out_path,"bim_overlap_analysis_results.RData",sep=""))

# Analyze the results
get_shared_loc_snpids<-function(x,y){
  inds = x[,5]==y[,5] & x[,6]==y[,6]
  return(cbind(x[inds,2],y[inds,2]))
}
final_shared_snps = cbind(shared_snps,shared_snps)
for(chr in names(intersected_locations)){
  final_shared_snps = rbind(final_shared_snps,get_shared_loc_snpids(
    intersected_locations[[chr]][[3]],intersected_locations[[chr]][[4]]
  ))
}
# Order the SNPs by location: prevents issues with PLINK
m = bim_data2[[1]][final_shared_snps[,2],]
ord = order(as.numeric(m[,1]),as.numeric(m[,4]))
final_shared_snps = final_shared_snps[ord,]
save(intersected_locations,final_shared_snps,file=paste(out_path,"bim_overlap_analysis_results.RData",sep=""))

load(paste(out_path,"bim_overlap_analysis_results.RData",sep=""))
extract_snps_using_plink(bfile2,final_shared_snps[,2],out_path,"file2_shared_snps","new_bed_2",
                        get_sh_default_prefix)
extract_snps_using_plink(bfile1,final_shared_snps[,1],out_path,"file1_shared_snps","new_bed_1",
                        get_sh_default_prefix)

# Print the new bed file for file 2
conv_ids = final_shared_snps[,1];names(conv_ids) = final_shared_snps[,2]
new_bed_2_bim = process_bim_data(paste(out_path,"new_bed_2",sep=""))
newbim2 = new_bed_2_bim[[1]]
newbim2[,2] = conv_ids[newbim2[,2]]
rownames(newbim2) = NULL
write.table(newbim2,file=paste(out_path,"new_bed_2_alt.bim",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)

# Convert to bgen
# August 2018: check rerunning with "id-delim=" because our ids have "_" in them
err_path = paste(out_path,"convert2bgen2.err",sep="")
log_path = paste(out_path,"convert2bgen2.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("~/apps/plink2/./plink2 --bed",paste(out_path,"new_bed_2.bed",sep=''),
                 "--bim",paste(out_path,"new_bed_2_alt.bim",sep=''),
                 "--fam", paste(out_path,"new_bed_2.fam",sep=''),
                 "--export bgen-1.3 --out",paste(out_path,"new_bgen_2",sep=''))
curr_sh_file = "convert2bgen2.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''), curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

err_path = paste(out_path,"convert2bgen1.err",sep="")
log_path = paste(out_path,"convert2bgen1.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("~/apps/plink2/./plink2 --bfile",paste(out_path,"new_bed_1",sep=''),
                 "--export bgen-1.3 --out",paste(out_path,"new_bgen_1",sep=''))
curr_sh_file = "convert2bgen1.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Run qctool
err_path = paste(out_path,"merge_qctool.err",sep="")
log_path = paste(out_path,"merge_qctool.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=32000,Ncpu=4,time="24:00:00")
curr_cmd = paste("~/apps/qctool_v2/build/release/./qctool_v2.0.1", 
            "-g",paste(out_path,"new_bgen_1.bgen",sep=''),
            "-s",paste(out_path,"new_bgen_1.sample",sep=''),
            "-g",paste(out_path,"new_bgen_2.bgen",sep=''),
            "-s",paste(out_path,"new_bgen_2.sample",sep=''),
            "-og",paste(out_path,"merged_data.bgen",sep=''),
            "-os",paste(out_path,"merged_data.sample",sep=''))
curr_sh_file = "merge_qctool_v2.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# convert merged data into bed
err_path = paste(out_path,"bgen_to_bed.err",sep="")
log_path = paste(out_path,"bgen_to_bed.log",sep="")
curr_cmd = paste("plink2 --bgen",paste(out_path,"merged_data.bgen",sep=""),
                 "--sample",paste(out_path,"merged_data.sample",sep=""),
                 "--make-bed --out",paste(out_path,"merged_data_qctool_bed",sep=''))
curr_sh_file = "bgen_to_bed.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000,plink_pkg = "plink/2.0a1"),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Add frq and PCA
err_path = paste(out_path,"merge_qctool_pca.err",sep="")
log_path = paste(out_path,"merge_qctool_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_qctool_bed",sep=''),
                 "--pca --freq --out",paste(out_path,"merged_data_qctool_bed",sep=''))
curr_sh_file = "merge_qctool_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# As a reference we do a simple merge using PLINK
err_path = paste(out_path,"merge_plink.err",sep="")
log_path = paste(out_path,"merge_plink.log",sep="")
curr_cmd = paste("plink --bed",paste(out_path,"new_bed_2.bed",sep=''),
                 "--bim",paste(out_path,"new_bed_2_alt.bim",sep=''),
                 "--fam", paste(out_path,"new_bed_2.fam",sep=''),
                 "--bmerge",paste(out_path,"new_bed_1",sep=''),
                 "--make-bed --out",paste(out_path,"merged_data_plink",sep=''))
curr_sh_file = "merge_plink.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Add frq and PCA
err_path = paste(out_path,"merge_plink_pca.err",sep="")
log_path = paste(out_path,"merge_plink_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_plink",sep=''),
                 "--pca 40 --freq --out",paste(out_path,"merged_data_plink",sep=''))
curr_sh_file = "merge_plink_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))






