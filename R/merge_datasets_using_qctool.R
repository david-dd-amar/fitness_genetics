# In this script we take two datasets of bed files (one per chromosome) and merge them
# using qctool v2
# Input datasets are preferably after running check_bim.pl
# For using qctool:
# 1. We use plink to transform the files to bgen
# 2. We merge the bgens using qctool
# 3. We convert the bgen to bed - the output
# 4. We flip snps to match the first dataset (-flip-to-match-cohort1)
# 5. We match snps by position and alleles (-compare-variants-by position,alleles)
# Output: merged bed and bgen files for each chromosome

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

dataset1 = "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/"
dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp_wo_jhu/impute2_1000gRef_out/check_bim_res/"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_1000g/"
system(paste("mkdir",out_path))

for(chr in 1:22){
  bfile1 = paste(dataset1,"chr",chr,sep="")
  bfile2 = paste(dataset2,"chr",chr,sep="")
  
  # # Some QC
  # bim1 = read.table(paste(bfile1,".bim",sep=""),stringsAsFactors = F)
  # bim2 = read.table(paste(bfile2,".bim",sep=""),stringsAsFactors = F)
  # rownames(bim1)=bim1[,2]
  # rownames(bim2)=bim2[,2]
  # which(table(bim1[,2])>1)
  # bim1[bim1[,2]=="rs10656307",]
  # length(intersect(bim1[,2],bim2[,2]))

  err_path = paste(out_path,"chr",chr,"_convert2bgen1.err",sep="")
  log_path = paste(out_path,"chr",chr,"_convert2bgen1.log",sep="")
  curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
    err_path,log_path,plink_pkg="plink/2.0a1",mem_size=32000,Ncpu=4)
  curr_cmd = paste("~/apps/plink2/./plink2 --bfile",bfile1,
                   "--export bgen-1.3 --out",
                   paste(out_path,"chr",chr,"_bgen_1",sep=''))
  curr_sh_file = paste(out_path,"chr",chr,"_convert2bgen1.sh",sep="")
  print_sh_file(curr_sh_file,curr_sbatch_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
  
  err_path = paste(out_path,"chr",chr,"_convert2bgen2.err",sep="")
  log_path = paste(out_path,"chr",chr,"_convert2bgen2.log",sep="")
  curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
    err_path,log_path,plink_pkg="plink/2.0a1",mem_size=32000,Ncpu=4)
  curr_cmd = paste("~/apps/plink2/./plink2 --bfile",bfile2,
                   "--export bgen-1.3 --out",
                   paste(out_path,"chr",chr,"_bgen_2",sep=''))
  curr_sh_file = paste(out_path,"chr",chr,"_convert2bgen2.sh",sep="")
  print_sh_file(curr_sh_file,curr_sbatch_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
  
  err_path = paste(out_path,"chr",chr,".err",sep="")
  log_path = paste(out_path,"chr",chr,".log",sep="")
  curr_sbatch_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(
    err_path,log_path,plink_pkg="plink/2.0a1",mem_size=64000,Ncpu=16,time="24:00:00")
  curr_cmd = paste("~/apps/qctool_v2/build/release/./qctool_v2.0.1", 
                   "-g",paste(out_path,"chr",chr,"_bgen_1.bgen",sep=''),
                   "-s",paste(out_path,"chr",chr,"_bgen_1.sample",sep=''),
                   "-g",paste(out_path,"chr",chr,"_bgen_2.bgen",sep=''),
                   "-s",paste(out_path,"chr",chr,"_bgen_2.sample",sep=''),
                   "-flip-to-match-cohort1 -compare-variants-by position,alleles",
                   "-threads 16",
                   "-og",paste(out_path,"chr",chr,"_merged_data.bgen",sep=''),
                   "-os",paste(out_path,"chr",chr,"_merged_data.sample",sep=''))
  curr_sh_file =  paste(out_path,"chr",chr,".sh",sep="")
  print_sh_file(curr_sh_file,curr_sbatch_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
  wait_for_job(120)
}