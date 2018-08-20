# In this script we take two bed files and merge them as follows:
# Input is two bed paths that are assumed to be comparable (no strand and MAF issues)
# 1. We use plink to transform the files to bgen
# 2. We merge the bgens using qctool
# 3. We convert the bgen to bed - the output
# 4. (Optional): given a case-control phenotype we check for issues using PLINK's flip scan
# As for the input, note that the analysis is not symmetric. Whenever possible, we convert the ids
# in file2 to match those in file1. For example if we use ukbb for file1 and Illumina ids for 
# file2 then non-standard ids such as (examXXX) will be mapped to ids from ukbb.

bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-updated"
bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/final_dataset_for_analysis-updated"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_with_ukbb1/"

check_bim_info = T
qctool_path = "/home/users/davidama/apps/qctool_v2/build/release/qctool_v2.0.1"

script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
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
save(intersected_locations,final_shared_snps,file=paste(out_path,"bim_overlap_analysis_results.RData",sep=""))

# Create the reduced files
reduce_snps_using_plink<-function(bfile,snps,out_path,snpfile,newbedfile,batch_script_func,...){
  # create the snps file
  write.table(t(t(as.character(snps))),
              file=paste(out_path,snpfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F)
  err_path = paste(out_path,"reduce_snps",snpfile,".err",sep="")
  log_path = paste(out_path,"reduce_snps",snpfile,".log",sep="")
  curr_cmd = paste("plink2 --bfile",bfile,
                   "--extract",paste(out_path,snpfile,".txt",sep=''),
                   "--freq --sort-vars --make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"reduce_snps",snpfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
load(paste(out_path,"bim_overlap_analysis_results.RData",sep=""))
reduce_snps_using_plink(bfile2,final_shared_snps[,2],out_path,"file2_shared_snps","new_bed_2",
                        get_sh_prefix_one_node_specify_cpu_and_mem,
                        plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
reduce_snps_using_plink(bfile1,final_shared_snps[,1],out_path,"file1_shared_snps","new_bed_1",
                        get_sh_prefix_one_node_specify_cpu_and_mem,
                        plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)

# Print the new bed file for file 2
conv_ids = final_shared_snps[,1];names(conv_ids) = final_shared_snps[,2]
new_bed_2_bim = process_bim_data(paste(out_path,new_bed_2,sep=""))
newbim2 = new_bed_2_bim[[1]]
newbim2[,2] = conv_ids[newbim2[,2]]
rownames(newbim2) = NULL
write.table(newbim2,file=paste(out_path,"new_bed_2_alt.bim",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)

# Convert to pgen
err_path = paste(out_path,"convert2bgen2.err",sep="")
log_path = paste(out_path,"convert2bgen2.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("plink2 --bed",paste(out_path,"new_bed_2.bed",sep=''),
                 "--bim",paste(out_path,"new_bed_2_alt.bim",sep=''),
                 "--fam", paste(out_path,"new_bed_2.fam",sep=''),
                 "--make-pgen --sort-vars --out",paste(out_path,"new_pgen2",sep=''))
curr_sh_file = "convert2bgen2.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''), curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

err_path = paste(out_path,"convert2bgen1.err",sep="")
log_path = paste(out_path,"convert2bgen1.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("plink2 --bfile",paste(out_path,"new_bed_1",sep=''),
                 "--make-pgen --sort-vars --out",paste(out_path,"new_pgen1",sep=''))
curr_sh_file = "convert2bgen1.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Convert to bgen
# August 2018: check rerunning with "id-delim=" because our ids have "_" in them
err_path = paste(out_path,"convert2bgen2.err",sep="")
log_path = paste(out_path,"convert2bgen2.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("~/apps/plink2/./plink2 --pfile",paste(out_path,"new_pgen2",sep=''),
                 "--export bgen-1.3 --out",paste(out_path,"new_bgen2",sep=''))
curr_sh_file = "convert2bgen2.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''), curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

err_path = paste(out_path,"convert2bgen1.err",sep="")
log_path = paste(out_path,"convert2bgen1.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("~/apps/plink2/./plink2 --pfile",paste(out_path,"new_pgen1",sep=''),
                 "--export bgen-1.3 --out",paste(out_path,"new_bgen1",sep=''))
curr_sh_file = "convert2bgen1.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Run qctool
err_path = paste(out_path,"merge_qctool.err",sep="")
log_path = paste(out_path,"merge_qctool.log",sep="")
curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
  err_path,log_path,plink_pkg="plink/2.0a1",mem_size=16000,Ncpu=1)
curr_cmd = paste("~/apps/qctool_v2/build/release/./qctool_v2.0.1", 
            "-g",paste(out_path,"new_bgen1.bgen",sep=''),
            "-s",paste(out_path,"new_bgen1.sample",sep=''),
            "-g",paste(out_path,"new_bgen2.bgen",sep=''),
            "-s",paste(out_path,"new_bgen2.sample",sep=''),
            "-og",paste(out_path,"merged_data.bgen",sep=''),
            "-os",paste(out_path,"merged_data.sample",sep=''))
curr_sh_file = "merge_qctool.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),curr_sbatch_prefix,curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))






