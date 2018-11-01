
# Define input for flow
script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/"
our_data_bed_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/1000g/merged_mega_data_autosomal"
map_files_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/1000g_data/1000GP_Phase3/"
shapeit_path = "/home/users/davidama/apps/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
impute2_path = "/home/users/davidama/apps/impute2/impute_v2.3.2_x86_64_static/impute2"
impute2_size = 5000000
ref_for_phasing = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp/1000g_data/1000GP_Phase3/"


# Split our data by chromosome
direct_geno_path = paste(job_dir,"direct_geno/",sep="")
# # Create a dir with bed file per chromosome
# system(paste("mkdir",direct_geno_path))
# setwd(direct_geno_path)
# for (j in 1:22){
#   err_path = paste("split",j,".err",sep="")
#   log_path = paste("split",j,".log",sep="")
#   curr_cmd = paste("plink --bfile",our_data_bed_path,
#                    "--chr",j,
#                    "--make-bed --out",paste(j,sep=''))
#   curr_sh_file = paste("split",j,".sh",sep="")
#   print_sh_file(curr_sh_file,get_sh_default_prefix(err_path,log_path),curr_cmd)
#   system(paste("sbatch",paste(curr_sh_file,sep='')))
# }

map_files = list.files(map_files_path)
map_files = map_files[grepl("^genetic_map",map_files)]
map_files = paste(map_files_path,map_files,sep="")
ref_files = list.files(ref_for_phasing)
hap_files = ref_files[grepl("hap.gz$",ref_files)]
hap_files = paste(ref_for_phasing,hap_files,sep="")
legend_files = ref_files[grepl("legend.gz$",ref_files)]
legend_files = paste(ref_for_phasing,legend_files,sep="")
sample_file = ref_files[grepl(".sample$",ref_files)]
sample_file = paste(ref_for_phasing,sample_file,sep="")

# # Run Shapeit for each chromosome, without a reference
# if(is.null(ref_for_phasing)){
#   shapeit_out_path = paste(job_dir,"shapeit_out/",sep="")
#   system(paste("mkdir",shapeit_out_path))
#   setwd(shapeit_out_path)
#   for (j in 1:22){
#     err_path = paste("shapeit_run_",j,".err",sep="")
#     log_path = paste("shapeit_run_",j,".log",sep="")
#     curr_bed = paste(direct_geno_path,j,".bed",sep="")
#     curr_bim = paste(direct_geno_path,j,".bim",sep="")
#     curr_fam = paste(direct_geno_path,j,".fam",sep="")
#     curr_gmap = map_files[grepl(paste("chr",j,"_",sep=""),map_files)]
#     curr_cmd = paste(shapeit_path,
#                      "--input-bed",curr_bed,curr_bim,curr_fam,
#                      "--input-map",curr_gmap,
#                      "-O",paste(j,"_phased",sep=''),
#                      "--thread 8 --seed 123456789")
#     curr_sh_file = paste("shapeit_run_",j,".sh",sep="")
#     print_sh_file(curr_sh_file,
#                   get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu = 8,mem_size = 16000),
#                   curr_cmd)
#     system(paste("sbatch",paste(curr_sh_file,sep='')))
#   }
# }

# Run Shapeit for each chromosome, with 1000G as a reference
if(!is.null(ref_for_phasing)){
  shapeit_1000gOut_path = paste(job_dir,"shapeit_1000gRef_out/",sep="")
  system(paste("mkdir",shapeit_1000gOut_path))
  setwd(shapeit_1000gOut_path)
  for (j in 1:22){
    curr_bed = paste(direct_geno_path,j,".bed",sep="")
    curr_bim = paste(direct_geno_path,j,".bim",sep="")
    curr_fam = paste(direct_geno_path,j,".fam",sep="")
    curr_hap = hap_files[grepl(paste("chr",j,"\\.",sep=""),hap_files)]
    curr_leg = legend_files[grepl(paste("chr",j,"\\.",sep=""),legend_files)]
    curr_gmap = map_files[grepl(paste("chr",j,"_",sep=""),map_files)]
    
    # # Step 3A from the tutorial: strand checks
    # err_path = paste("shapeit_run_check_",j,".err",sep="")
    # log_path = paste("shapeit_run_check_",j,".log",sep="")
    # curr_cmd = paste(shapeit_path,
    #                  "-check --input-bed",curr_bed,curr_bim,curr_fam,
    #                  "--input-map",curr_gmap,
    #                  "--input-ref", curr_hap,curr_leg,sample_file,
    #                  "--output-log",paste(j,".alignments",sep=''),
    #                  "--thread 8 --seed 123456789")
    # curr_sh_file = paste("shapeit_run_check_",j,".sh",sep="")
    # print_sh_file(curr_sh_file,
    #               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu = 8,mem_size = 16000),
    #               curr_cmd)
    # system(paste("sbatch",paste(curr_sh_file,sep='')))
    # wait_for_job(waittime = 120)
    
    # Step 3B: check if strand issues where found
    strand_files = list.files(shapeit_1000gOut_path)
    strand_files = strand_files[grepl(paste(j,".alignments.strand",sep=''),strand_files)]
    exclude_line = ""
    if(length(strand_files)>0){
      strand_files = strand_files[grepl(paste(j,".alignments.strand.exclude",sep=''),strand_files)]
      exclude_line = paste("--exclude-snp",strand_files)
    }
    
    # Final step: run Shapeit
    err_path = paste("shapeit_run_",j,".err",sep="")
    log_path = paste("shapeit_run_",j,".log",sep="")
    curr_cmd = paste(shapeit_path,
                     "--input-bed",curr_bed,curr_bim,curr_fam,
                     "--input-map",curr_gmap,exclude_line,
                     "--input-ref", curr_hap,curr_leg,sample_file,
                     "-O",paste(j,"_phased",sep=''),
                     "--thread 8 --seed 123456789")
    curr_sh_file = paste("shapeit_run_",j,".sh",sep="")
    print_sh_file(curr_sh_file,
                  get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu = 8,mem_size = 16000),
                  curr_cmd)
    system(paste("sbatch",paste(curr_sh_file,sep='')))
  }
}

# Option 2: shapeit with 1000G as ref
shapeit_out_path = paste(job_dir,"shapeit_1000gRef_out/",sep="")
impute2_out_path = paste(job_dir,"impute2_1000gRef_out/",sep="")
system(paste("mkdir",impute2_out_path))
setwd(impute2_out_path)
hap_files = ref_files[grepl("hap.gz$",ref_files)]
hap_files = paste(ref_for_phasing,hap_files,sep="")
legend_files = ref_files[grepl("legend.gz$",ref_files)]
legend_files = paste(ref_for_phasing,legend_files,sep="")
tmp_chunck_files_per_chr = list()
for (j in 1:22){
  curr_haps = paste(shapeit_out_path,j,"_phased.haps",sep='')
  curr_gmap = map_files[grepl(paste("chr",j,"_",sep=""),map_files)]
  curr_ref_hap = hap_files[grepl(paste("chr",j,"\\.",sep=""),hap_files)]
  curr_leg = legend_files[grepl(paste("chr",j,"\\.",sep=""),legend_files)]
  curr_gmap_info = read.table(curr_gmap,stringsAsFactors = F)
  max_chr_ind = as.numeric(curr_gmap_info[nrow(curr_gmap_info),1])
  start_ind = 1
  chunk_count = 1
  chunk2out = c()
  while(start_ind <= max_chr_ind){
    err_path = paste("tmp_chr",j,"_chunck",chunk_count,".err",sep="")
    log_path = paste("tmp_chr",j,"_chunck",chunk_count,".log",sep="")
    curr_out = paste("chr",j,"_chunck",chunk_count,".imputed",sep="")
    curr_out_i = paste("chr",j,"_chunck",chunk_count,"_info",sep="")
    chunk2out = c(chunk2out,curr_out,curr_out_i)
    if(is.element(curr_out,set=list.files())){
      print(paste("Skipping:",curr_out))
      chunk_count = chunk_count + 1
      start_ind = start_ind + impute2_size
      next
    }
    # create all jobs for the current chromosome
    curr_cmd = paste(impute2_path, "-use_prephased_g",
                     "-known_haps_g",curr_haps,
                     "-h",curr_ref_hap,
                     "-l",curr_leg,
                     "-m",curr_gmap,
                     "-int",start_ind,(start_ind + impute2_size -1),
                     "-Ne 20000",
                     "-o", curr_out,
                     "-i",curr_out_i
    )
    curr_sh_file = paste("tmp_chr",j,"_chunck",chunk_count,".sh",sep="")
    print_sh_file(curr_sh_file,
                  get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu = 4,mem_size = 32000),
                  curr_cmd)
    system(paste("sbatch",paste(curr_sh_file,sep='')))
    chunk_count = chunk_count + 1
    start_ind = start_ind + impute2_size
  }
  # For each chr: wait for all imp jobs to finish
  tmp_chunck_files_per_chr[[j]] = chunk2out
  save(tmp_chunck_files_per_chr,file="tmp_chunck_files_per_chr.RData")
}
wait_for_job(waittime = 60)
wait_for_job(waittime = 60)

## Organazie the output

# 1. Put all info files in a dir
info_dir = "info_files/"
system(paste("mkdir",info_dir))
allfiles = list.files()
info_files = allfiles[grepl("info",allfiles)]
for(ff in info_files){system(paste("mv",ff,info_dir))}
info_files = paste(info_dir,info_files,sep="")
info_files = info_files[!grepl("by_sample",info_files)]
info_file_summary = paste(info_dir,"our_data_info.txt",sep="")
for(i in 1:length(info_files)){
  cmd = paste("less",info_files[i],"| grep --perl \"^\\d\" > ",info_file_summary)
  if(i>1){
    cmd = paste("less",info_files[i],"| grep --perl \"^\\d\" >> ",info_file_summary)
  }
  system(cmd)
}
info_example = read.table(info_files[1],header=T)
our_data_info = read.table(info_file_summary)
colnames(our_data_info) = colnames(info_example)
save(our_data_info,file="our_data_info.RData")

# 2. Concatenate all impute result files
load("tmp_chunck_files_per_chr.RData")
for(j in 1:22){
  print(j)
  chunk2out = tmp_chunck_files_per_chr[[j]]
  chunk2out = chunk2out[grepl("imputed$",chunk2out)]
  chunk2out = intersect(chunk2out,list.files(getwd()))
  if(length(chunk2out)==0){next}
  # put all imputed in one file
  system(paste(paste(c("cat",chunk2out),collapse = " "),paste("> chr",j,".imputed",sep="")))
  for(chunk in chunk2out){
    system(paste("rm",chunk))
  }
}

# 3. Put all other files in archive dirs
for(j in 1:22){
  print(j)
  # Clean the directory
  run_archive_dir = paste("chr",j,"_run_archive",sep="")
  system(paste("mkdir",run_archive_dir))
  rr = list.files(getwd())
  rr = rr[grepl("chunck",rr) & grepl(paste("chr",j,"_",sep=""),rr)]
  for(ff in rr){system(paste("mv",ff,run_archive_dir))}
}

# 4. gzip the imputation files
for(j in 1:22){
  print(j)
  system(paste("gzip",paste("chr",j,".imputed &",sep="")))
}

# 5. create the sample file
fam = read.table(paste(our_data_bed_path,".fam",sep=""))
sample_info = as.matrix(fam[,c(1,2)])
sample_info = cbind(sample_info,rep(0,nrow(sample_info)))
sample_info = cbind(sample_info,fam[,5])
l1 = c("ID_1","ID_2","missing","sex")
l2 = c("0","0","0","D")
sample_info = rbind(l1,l2,sample_info)
write.table(sample_info,col.names = F,row.names = F,quote = F,sep=" ",file="sample_file.sample")

# 6. Transform to bed files
for (j in 1:22){
  err_path = paste("impute2bed_",j,".err",sep="")
  log_path = paste("impute2bed_",j,".log",sep="")
  curr_cmd = paste("plink --gen",paste("chr",j,".imputed.gz",sep=""),
                   "--sample sample_file.sample",
                   "--oxford-single-chr",j,
                   "--make-bed --out",paste("chr",j,sep=""))
  curr_sh_file = paste("impute2bed_",j,".sh",sep="")
  print_sh_file(curr_sh_file,get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",curr_sh_file))
}

# 7. Go over the info and qc of the imputation
# Exclude SNPs with low MAF
for (j in 1:22){
  err_path = paste("maf_",j,".err",sep="")
  log_path = paste("maf_",j,".log",sep="")
  curr_cmd = paste("plink --bfile",paste("chr",j,sep=""),
                   "--maf 0.001 --freq",
                   "--make-bed --out",paste("chr",j,sep=""))
  curr_sh_file = paste("maf_",j,".sh",sep="")
  print_sh_file(curr_sh_file,get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
wait_for_job(waittime = 120)

# Exclude SNPs that are from our data and have low concodrance score
# Exclude imputed SNPs with a low certainty score
load(paste(info_dir,"our_data_info.RData",sep=""))
rownames(our_data_info) = our_data_info$rs_id
# # Examine our concordance scores for "zero" pval snps
# gwas_res = read.table(
#           "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/gwas/ukbb_vs_cooper_logistic.ExerciseGroup.glm.logistic.hybrid.adjusted",
#           stringsAsFactors = F,header=F)
# rownames(gwas_res) = gwas_res[,2]
# inds = intersect(rownames(gwas_res),rownames(our_data_info) )
# x1 = -log(gwas_res[inds,3],base=10)
# x2 = our_data_info[inds,]$concord_type0
our_low_concord_snps = as.character(our_data_info$rs_id[our_data_info$concord_type0 < 0.8])
info_files = list.files(info_dir)
info_files = paste(info_dir,info_files,sep="")
info_files = info_files[!grepl("by_sample",info_files)]
low_cert_snps = c()
for(ff in info_files){
  curr_data = read.table(ff,header=T,stringsAsFactors = F)
  curr_excluded = curr_data$rs_id[curr_data$certainty < 0.9]
  low_cert_snps = c(low_cert_snps,curr_excluded)
  print(length(low_cert_snps))
}
save(our_data_info,low_cert_snps,our_low_concord_snps,file=paste(info_dir,"our_data_info.RData",sep=""))

to_rem = c(low_cert_snps,our_low_concord_snps)
to_rem_snp_file = "impute_qc_snps_to_remove"
for(j in 1:22){
  bfile = paste("chr",j,sep="")
  curr_out_path = paste(getwd(),'/',sep="")
  curr_analysis_name = paste(to_rem_snp_file,j,sep="")
  exclude_snps_using_plink(bfile,to_rem,curr_out_path,curr_analysis_name,bfile)
}

# 8. Run check bim vs 1000G before merge
# 8.1 update the frq files
for (j in 1:22){
  err_path = paste("maf_",j,".err",sep="")
  log_path = paste("maf_",j,".log",sep="")
  curr_cmd = paste("plink --bfile",paste("chr",j,sep=""),
                   "--freq",
                   "--out",paste("chr",j,sep=""))
  curr_sh_file = paste("maf_",j,".sh",sep="")
  print_sh_file(curr_sh_file,get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
wait_for_job(waittime = 120)

# 8.2 Run check bim
for (j in 1:22){
  setwd(impute2_out_path)
  curr_dir = paste(impute2_out_path,"check_bim_chr",j,'/',sep="")
  bedfile = paste(impute2_out_path,"chr",j,".bim",sep="")
  freqfile = paste(impute2_out_path,"chr",j,".frq",sep="")
  run_check_bim_analysis(curr_dir,bedfile,freqfile)
}
wait_for_job(waittime = 120)
for (j in 1:22){
  setwd(impute2_out_path)
  curr_dir = paste(impute2_out_path,"check_bim_chr",j,'/',sep="")
  bedfile = paste(impute2_out_path,"chr",j,sep="")
  bedfile_short = paste("chr",j,sep="")
  run_check_bim_output_script(curr_dir,bedfile_short,bedfile)
}
wait_for_job(waittime = 120)

