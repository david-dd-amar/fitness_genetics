
# Input
# imputed
# external_files_path = "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/" # imputed
# analysis_name = "ukbb_imputed_20k_rand_controls_sex_age"
# direct
external_files_path = "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/" # direct
analysis_name = "ukbb_direct_20k_rand_controls_sex_age"
# metadata
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age.txt"
external_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/20k_rand_controls_sex_age_with_info.txt"

# set the output dir
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/"
chrs = 1:22
all_files = list.files(external_files_path)
if(analysis_name != ""){
  out_path = paste(out_path,analysis_name,"/",sep="")
  system(paste("mkdir",out_path))
}

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

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
wait_for_job(waittime = 120)

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
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(120)

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
# Transform the dataset into HRC-based data
setwd(out_path)
# Run the check_bim analysis
err_path = paste(out_path,"run_check_bim.err",sep="")
log_path = paste(out_path,"run_check_bim.log",sep="")
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",out_path))
curr_cmd = paste("perl", paste(out_path, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(out_path,"merged_control_geno.bim",sep=''),
                 "-f", paste(out_path,"merged_control_geno.frq",sep=''),
                 "-hrc -p EU -t 0.1 -r",
                 "/home/users/davidama/apps/check_bim/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
curr_sh_file = "run_check_bim.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
system(paste("less ",out_path,"Run-plink.sh | grep TEMP > ",out_path,"Run-plink_hrc.sh",sep=""))
run_sh_lines = readLines(paste(out_path,"Run-plink_hrc.sh",sep=""))
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "-hrc_updated")
write.table(file=paste(out_path,"Run-plink_hrc.sh",sep=""),t(t(run_sh_lines)),
            quote=F,row.names = F,col.names = F)
# Run check-bim's output script
err_path = paste(out_path,"run_check_bim_update.err",sep="")
log_path = paste(out_path,"run_check_bim_update.log",sep="")
plink_commands = readLines(paste(out_path,"Run-plink_hrc.sh",sep=""))
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),plink_commands)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# Use 1000G instead of HRC
setwd(out_path)
# Run the check_bim analysis
err_path = paste(out_path,"run_check_bim2.err",sep="")
log_path = paste(out_path,"run_check_bim2.log",sep="")
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",out_path))
curr_cmd = paste("perl", paste(out_path, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(out_path,"merged_control_geno.bim",sep=''),
                 "-f", paste(out_path,"merged_control_geno.frq",sep=''),
                 "-1000g -p EUR -t 0.1 -r ",
                 "/home/users/davidama/apps/check_bim/1000GP_Phase3_combined.legend")
curr_sh_file = "run_check_bim2.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=64000),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(120)
system(paste("less ",out_path,"Run-plink.sh | grep TEMP > ",out_path,"Run-plink_1000g.sh",sep=""))
run_sh_lines = readLines(paste(out_path,"Run-plink_1000g.sh",sep=""))
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "-1000g_updated")
write.table(file=paste(out_path,"Run-plink_1000g.sh",sep=""),t(t(run_sh_lines)),
            quote=F,row.names = F,col.names = F)
# Run check-bim's output script
err_path = paste(out_path,"run_check_bim_update.err",sep="")
log_path = paste(out_path,"run_check_bim_update.log",sep="")
plink_commands = readLines(paste(out_path,"Run-plink_1000g.sh",sep=""))
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),plink_commands)
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=64000),plink_commands)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))

# ## MAF filter to reduce file size
# err_path = paste(out_path,"merge_control_beds.err",sep="")
# log_path = paste(out_path,"merge_control_beds.log",sep="")
# curr_cmd = paste("plink --bfile",all_out_bed_files[1],
#                  "--merge-list",allfiles_path,
#                  "--make-bed --out",paste(out_path,"merged_control_geno",sep=''))
# curr_sh_file = "merged_control_beds.sh"
# print_sh_file(paste(out_path,curr_sh_file,sep=''),
#               get_sh_prefix_bigmem(err_path,log_path,Ncpu=1,mem_size=256000),curr_cmd)
# system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))




