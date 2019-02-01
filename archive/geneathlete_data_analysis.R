# Analysis of the GENEATHLETE dir
WD = "/oak/stanford/groups/euan/projects/fitness_genetics/GENEATHLETE/"
setwd(WD)
job_dir=WD
all_files = paste(getwd(),'/',list.files(),sep='')
# load functions
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Check _1 vs _2 : conclusion - just split, can merge
# map_1_files = sort(all_files[grepl("_1.map",all_files)])
# map_2_files = sort(all_files[grepl("_2.map",all_files)])
# for(j in 1:length(map_1_files)){
#   f1 = map_1_files[j]
#   f2 = map_2_files[j]
#   if(gsub("_\\d.map","",f1)!=gsub("_\\d.map","",f2)){print(paste("ERROR in file number:",j,f1,f2))}
#   d1 = read.table(f1,header = F)
#   d2 = read.table(f2,header = F)
#   print(length(setdiff(d1$V2,d2$V2)))
#   print(c(nrow(d1),nrow(d2)))
# }

####################################################################################################
####################################################################################################
####################################################################################################
# merge all geneathlete data into one file
map_files = sort(all_files[grepl(".map$",all_files)])
m = c()
for(j in 1:length(map_files)){
  f_map = map_files[j]
  f_ped = gsub("map","ped",f_map)
  m = rbind(m,c(f_ped,f_map))
}
write.table(m[-1,],file="allfiles.txt",sep=" ",quote=F,row.names = F,col.names = F)

jobs_before = get_my_jobs()
err_path = paste(job_dir,"merged_peds.err",sep="")
log_path = paste(job_dir,"merged_peds.log",sep="")
first_re = gsub(".ped","",m[1,1])
curr_cmd = paste("plink --file",first_re,
                 "--merge-list",paste(job_dir,"allfiles.txt",sep=""),
                 "--missing --recode --out",paste(job_dir,"merged_peds",sep=''))
curr_sh_file = "merge_preds.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(jobs_before,5)
list.files(job_dir)




