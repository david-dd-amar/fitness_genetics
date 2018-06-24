
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

external_files_path = "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_dir_genotype/"
external_control_ids = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/10k_rand_controls_sex_age.txt"
our_bed_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/maf_filter"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/"
chrs = 1:22
all_files = list.files(external_files_path)
analysis_name = "ukbb_10k_rand_controls_sex_age"

####################################################################################################
####################################################################################################
####################################################################################################
# Create control files
for (chr in chrs){
  curr_file = all_files[grepl(".bed$",all_files) & grepl(paste("chr",chr,"_",sep=""),all_files)]
  curr_file = gsub(pattern = ".bed$",replacement = "",curr_file)
  err_path = paste(out_path,analysis_name,".err",sep="")
  log_path = paste(out_path,analysis_name,".log",sep="")
  curr_cmd = paste("plink --bfile",paste(external_files_path,curr_file,sep=''),
                   "--keep",external_control_ids,
                   "--make-bed --out",paste(out_path,analysis_name,"_chr",chr,sep=""))
  curr_sh_file = "get_controls.sh"
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}


####################################################################################################
####################################################################################################
####################################################################################################
# Merge bed/ped files

all_out_bed_files = list.files(out_path)
all_out_bed_files = all_out_bed_files[grepl(analysis_name,all_out_bed_files) &
                                        grepl(".bed$",all_out_bed_files)]
all_out_bed_files = gsub(".bed","",all_out_bed_files)
allfiles_path = paste(out_path,analysis_name,"_allfiles.txt",sep="")
write.table(t(t(all_out_bed_files[-1])),
            file = allfiles_path,sep="",row.names = F,col.names = F,quote = F)

jobs_before = get_my_jobs()
err_path = paste(out_path,"merge_control_beds.err",sep="")
log_path = paste(out_path,"merge_control_beds.log",sep="")
curr_cmd = paste("plink --bfile",paste(out_path,all_out_bed_files[1],sep=''),
                 "--merge-list",allfiles_path,
                 "--make-bed --out",paste(out_path,analysis_name,"merged_control_geno",sep=''))
curr_sh_file = "merged_control_beds.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
list.files(job_dir)

all_out_bed_files = list.files(out_path)
all_out_bed_files = all_out_bed_files[grepl(analysis_name,all_out_bed_files) &
                                        grepl(".bed$",all_out_bed_files)]
readLines(log_path)

# merge with our bed file
controls_bed = paste(out_path,analysis_name,"merged_control_geno",sep='')
jobs_before = get_my_jobs()
err_path = paste(out_path,"merge_with_our_bed.err",sep="")
log_path = paste(out_path,"merge_with_our_bed.log",sep="")
curr_cmd = paste("plink --bfile",our_bed_path,
                 "--bmerge",controls_bed,
                 "--make-bed --out",paste(out_path,analysis_name,"_merged_bed_final",sep=''))
curr_sh_file = "merge_with_our_bed.sh"
print_sh_file(paste(out_path,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
list.files(job_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# Create phe file

####################################################################################################
####################################################################################################
####################################################################################################
# Run GWAS
