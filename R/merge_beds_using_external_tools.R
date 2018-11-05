# In this script we take two bed files and merge them as follows:
# Input is two bed paths that are assumed to be comparable (no strand and MAF issues)
# For using qctool:
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

# # September 2018 1: new MEGA analysis, HRC as the panel
# bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-hrc_updated"
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data_autosomal-hrc_updated"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_hrc/"
# 
# # September 2018 2: new MEGA analysis, 1000 genomes as the panel
# bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated"
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data_autosomal-updated"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g/"
# 
# # September 2018 2: new MEGA analysis with PCA filter, 1000 genomes as the panel, sanity check without JHU 
# bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated"
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/1000g/merged_mega_data_autosomal_after_maf_after_pca"
# out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/"

# Nov 2018: Merge imputed data
bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated"
bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/impute2_1000gRef_out/check_bim_res/merged_geno"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/with_ukbb/"

try(system(paste("mkdir",out_path)))

check_bim_info = T
qctool_path = "/home/users/davidama/apps/qctool_v2/build/release/qctool_v2.0.1"

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

remove_JHU = grepl("sanity",out_path)
maf_for_pca = 0.05

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# maf1 = 0.01
# maf2 = 0.01 # our data was already filtered in our prepro script
# 
# # create maf reduced copies of the data
# err_path = paste(out_path,"maf_filter1.err",sep="")
# log_path = paste(out_path,"maf_filter1.log",sep="")
# curr_cmd = paste("plink --bfile",bfile1,
#                  "--maf",maf1,
#                  "--make-bed --freq",
#                  "--out",paste(out_path,"bfile1",sep=""))
# curr_sh_file = "maf_filter1.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# 
# err_path = paste(out_path,"maf_filter2.err",sep="")
# log_path = paste(out_path,"maf_filter2.log",sep="")
# curr_cmd = paste("plink --bfile",bfile2,
#                  "--maf",maf2,
#                  "--make-bed --freq",
#                  "--out",paste(out_path,"bfile2",sep=""))
# curr_sh_file = "maf_filter2.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# wait_for_job()
# 
# bfile1 = paste(out_path,"bfile1",sep="")
# bfile2 = paste(out_path,"bfile2",sep="")
# 
####################################################################################################
####################################################################################################
####################################################################################################
# Compare the bim files
# (1) Check SNP intersect, locations, and which snps must be flipped before we analyze
bim_data1 = read.table(paste(bfile1,".bim",sep=""),stringsAsFactors = F,header = F)
bim_data2 = read.table(paste(bfile2,".bim",sep=""),stringsAsFactors = F,header = F)
rownames(bim_data2) = bim_data2[,2]
rownames(bim_data1) = bim_data1[,2]

if(remove_JHU){
  print(paste("num variants before JHU removal:",nrow(bim_data2[[1]])))
  bim_data2[[1]] = bim_data2[[1]][!grepl("JHU",bim_data2[[1]][,2]),]
  print(paste("num variants after JHU removal:",nrow(bim_data2[[1]])))
}

shared_snps = intersect(bim_data2[,2],bim_data1[,2])
length(shared_snps)
intersected_locations = list()
for(chr in unique(bim_data2[,1])){
  print(paste("chromosome:",chr))
  rows1 = bim_data1[,1]==chr
  rows2 = bim_data2[,1]==chr
  x1 = bim_data1[rows1,]
  x2 = bim_data2[rows2,]
  curr_loc_intersect = intersect(x1[,4],x2[,4])
  alt_ids1 = apply(x1[,4:6],1,paste,collapse=";")
  alt_ids2 = apply(x2[,4:6],1,paste,collapse=";")
  alt_ids2_sanity = apply(x2[,c(4,6,5)],1,paste,collapse=";")
  print(paste("intersection:",length(intersect(alt_ids1,alt_ids2))))
  print(paste("number of shared locations:",length(curr_loc_intersect)))
  print(paste("intersection sanity check (should be zero):",length(intersect(alt_ids1,alt_ids2_sanity))))
  alt_intersect = intersect(alt_ids1,alt_ids2)
  rownames(x1) = alt_ids1
  rownames(x2) = alt_ids2
  x1 = x1[alt_intersect,];x2 = x2[alt_intersect,]
  intersected_locations[[chr]] = list(x1,x2)
}
save(intersected_locations,file=paste(out_path,"bim_overlap_analysis_results.RData",sep=""))

# Analyze the results
final_shared_snps = c()
for(chr in 1:length(intersected_locations)){
  final_shared_snps = rbind(final_shared_snps,
    cbind(intersected_locations[[chr]][[1]][,2],intersected_locations[[chr]][[2]][,2])
  )
}
# Order the SNPs by location: prevents issues with PLINK
m = bim_data2[final_shared_snps[,2],]
ord = order(as.numeric(m[,1]),as.numeric(m[,4]))
final_shared_snps = final_shared_snps[ord,]
save(intersected_locations,final_shared_snps,file=paste(out_path,"bim_overlap_analysis_results.RData",sep=""))

load(paste(out_path,"bim_overlap_analysis_results.RData",sep=""))
extract_snps_using_plink(bfile2,final_shared_snps[,2],out_path,"file2_shared_snps","new_bed_2",
                        get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
extract_snps_using_plink(bfile1,final_shared_snps[,1],out_path,"file1_shared_snps","new_bed_1",
                        get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
wait_for_job(waittime = 60)

# Compare the new two beds
new_bed_2_bim = read.table(paste(out_path,"new_bed_2.bim",sep=""))
rownames(new_bed_2_bim) = new_bed_2_bim[,2]
new_bed_1_bim = read.table(paste(out_path,"new_bed_1.bim",sep=""))
rownames(new_bed_1_bim) = new_bed_1_bim[,2]
table(new_bed_1_bim[,1]==new_bed_2_bim[,1])
table(new_bed_1_bim[,4]==new_bed_2_bim[,4])
table(as.character(new_bed_1_bim[,5])==as.character(new_bed_2_bim[,5]))
table(as.character(new_bed_1_bim[,6])==as.character(new_bed_2_bim[,6]))
diff_rows = which(as.character(new_bed_1_bim[,6])!=as.character(new_bed_2_bim[,6]))
new_bed_1_bim[diff_rows[1:10],]
new_bed_2_bim[diff_rows[1:10],]

# Print the new bed file for file 2
conv_ids = final_shared_snps[,1];names(conv_ids) = final_shared_snps[,2]
new_bed_2_bim[,2] = conv_ids[new_bed_2_bim[,2]]
rownames(new_bed_2_bim) = NULL
write.table(new_bed_2_bim,file=paste(out_path,"new_bed_2_alt.bim",sep=""),sep="\t",
            row.names=F,col.names=F,quote=F)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# Merge using qctool
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
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              Ncpu=4,mem_size=64000,plink_pkg = "plink/2.0a1"),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Add frq and PCA
# Snp prune (with MAF)
analysis_name = "qctool_prune"
err_path = paste(out_path,analysis_name,"_ld_report.err",sep="")
log_path = paste(out_path,analysis_name,"_ld_report.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_qctool_bed",sep=''),
                 "--indep-pairwise 250 10",0.1,
                 "--out",paste(out_path,analysis_name,sep=""))
curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()
# Run PCA
err_path = paste(out_path,"merge_qctool_pca.err",sep="")
log_path = paste(out_path,"merge_qctool_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_qctool_bed",sep=''),
                 "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
                 "--pca 40 --out",paste(out_path,"merged_data_qctool_bed",sep=''))
curr_sh_file = "merge_qctool_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# Run relatedness
err_path = paste(out_path,"merge_qctool_relatedness.err",sep="")
log_path = paste(out_path,"merge_qctool_relatedness.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_qctool_bed",sep=''),
                 "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
                 "--genome --min 0.2 --out",paste(out_path,"merged_data_qctool_bed",sep=''))
curr_sh_file = "merge_qctool_relatedness.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# Add freq on all SNPs
err_path = paste(out_path,"merge_qctool_frq.err",sep="")
log_path = paste(out_path,"merge_qctool_frq.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_qctool_bed",sep=''),
                 "--freq --out",paste(out_path,"merged_data_qctool_bed",sep=''))
curr_sh_file = "merge_qctool_frq.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

###############################
# A simple merge using PLINK, the results should be similar to qctool
###############################
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

# Snp prune
analysis_name = "plink_prune"
err_path = paste(out_path,analysis_name,"_ld_report.err",sep="")
log_path = paste(out_path,analysis_name,"_ld_report.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_plink",sep=''),
                 "--indep-pairwise 250 10",0.1,
                 "--out",paste(out_path,analysis_name,sep=""))
curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job()
print(paste("number of snps after prunning:",
            length(readLines(paste(out_path,analysis_name,".prune.in",sep="")))))
# Run PCA
err_path = paste(out_path,"merge_plink_pca.err",sep="")
log_path = paste(out_path,"merge_plink_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_plink",sep=''),
                 "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
                 "--pca 40 --out",paste(out_path,"merged_data_plink",sep=''))
curr_sh_file = "merge_plink_pca.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# Run relatedness
err_path = paste(out_path,"merge_plink_relatedness.err",sep="")
log_path = paste(out_path,"merge_plink_relatedness.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_plink",sep=''),
                 "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
                 "--genome --min 0.2 --out",paste(out_path,"merged_data_plink",sep=''))
curr_sh_file = "merge_plink_relatedness.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
# Freq on all snps
err_path = paste(out_path,"merge_plink_frq.err",sep="")
log_path = paste(out_path,"merge_plink_frq.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_data_plink",sep=''),
                 "--freq --out",paste(out_path,"merged_data_plink",sep=''))
curr_sh_file = "merge_plink_frq.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# run_prune_pca_rl_analysis<-function(bfile,analysis_name,out_path,maf=0.05,){
#   # Snp prune
#   err_path = paste(out_path,analysis_name,"_ld_report.err",sep="")
#   log_path = paste(out_path,analysis_name,"_ld_report.log",sep="")
#   curr_cmd = paste("plink --bfile",bfile,
#                    "--indep-pairwise 250 10",0.1,
#                    "--out",paste(out_path,analysis_name,sep=""))
#   curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
#   print_sh_file(paste(out_path,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
#   system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#   wait_for_job()
#   # Run PCA
#   err_path = paste(out_path,"merge_plink_pca.err",sep="")
#   log_path = paste(out_path,"merge_plink_pca.log",sep="")
#   curr_cmd = paste("plink --bfile",bfile,
#                    "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
#                    "--pca 40 --out",bfile)
#   curr_sh_file = "merge_plink_pca.sh"
#   print_sh_file(paste(out_path,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
#   system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#   # Run relatedness
#   err_path = paste(out_path,"merge_plink_relatedness.err",sep="")
#   log_path = paste(out_path,"merge_plink_relatedness.log",sep="")
#   curr_cmd = paste("plink --bfile",bfile,
#                    "--extract", paste(out_path,analysis_name,".prune.in",sep=""),
#                    "--genome --min 0.2 --out",bfile)
#   curr_sh_file = "merge_plink_relatedness.sh"
#   print_sh_file(paste(out_path,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
#   system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#   # Freq on all snps
#   err_path = paste(out_path,"merge_plink_frq.err",sep="")
#   log_path = paste(out_path,"merge_plink_frq.log",sep="")
#   curr_cmd = paste("plink --bfile",bfile,
#                    "--freq --out",bfile)
#   curr_sh_file = "merge_plink_frq.sh"
#   print_sh_file(paste(out_path,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
#   system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
#   
# }




