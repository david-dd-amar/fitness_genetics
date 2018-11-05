####################################################################################################
####################################################################################################
####################################################################################################
######## Some useful functions
get_sh_default_prefix<-function(err="",log=""){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=6:00:00",
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=2",
      "#SBATCH --mem=16000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}

# plink2: plink/2.0a1
get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",plink_pkg = "plink/1.90b5.3",Ncpu,mem_size,time="6:00:00"){
  partition_line = "#SBATCH --partition=euan,mrivas,normal,owners"
  if(mem_size>128000){
    partition_line = "#SBATCH --partition=bigmem,euan,mrivas"
  }
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=",time,sep=""),
      partition_line,
      "#SBATCH --nodes=1",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      paste("module load",plink_pkg)
    )
  )
}

get_sh_prefix_bigmem<-function(err="",log="",plink_pkg = "plink/1.90b5.3",Ncpu=1,mem_size,time="6:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=",time,sep=""),
      "#SBATCH --partition=bigmem",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      paste("module load",plink_pkg)
    )
  )
}

print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}

get_my_jobs<-function(){
  tmp = paste("tmp",abs(rnorm(1)),sep='')
  system(paste("sacct > ",tmp),wait = T)
  jobs = readLines(tmp)[-c(1:2)]
  jobs = jobs[!grepl("^\\d+\\.",jobs)]
  jobs = t(sapply(jobs,function(x)strsplit(x,split="\\s+")[[1]]))
  rownames(jobs) = as.character(jobs[,1])
  jobs = jobs[,-1]
  jobs = jobs[jobs[,1]!="bash",]
  system(paste("rm",tmp))
  new_jobs = rownames(jobs)[jobs[,5]=="RUNNING" | jobs[,5]=="PENDING"]
  if(length(new_jobs)==1){
    return(matrix(jobs[new_jobs,],nrow=1))
  }
  return(jobs[new_jobs,])
}
get_job_id<-function(x){return(x[1])}
wait_for_job<-function(jobs_before=NULL,waittime=30,max_wait=6000){
  Sys.sleep(waittime)
  curr_jobs = get_my_jobs()
  if(length(curr_jobs)==0){return(NULL)}
  new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
  if(length(new_jobs)==0){return(matrix(NA,0,0))}
  print(paste("new added jobs are: ",new_jobs))
  i = 1
  while(T){
    Sys.sleep(waittime)
    curr_jobs = get_my_jobs()
    if(length(curr_jobs)==0){return(NULL)}
    new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
    if(length(new_jobs)==0){return(matrix(NA,0,0))}
    i=i+1
    if(i>max_wait){break}
  }
}

correct_dups_in_sample_metadata<-function(x){
  nns = apply(x[,1:2],1,paste,collapse="_")
  dups = names(which(table(nns)>1))
  to_keep = rep(F,length(nns))
  for(i in 1:length(nns)){
    if(!is.element(nns[i],set=dups)){to_keep[i]=T;next}
    curr_inds = which(nns==nns[i])
    num_nas = apply(is.na(x[curr_inds,]),1,sum)
    curr_inds = curr_inds[num_nas==min(num_nas)]
    to_keep[curr_inds[1]]=T
  }
  x = x[to_keep,]
  rownames(x) = nns[to_keep]
  return(x)
}
read_plink_table<-function(path,has_header=T,...){
  y = read.delim(path,stringsAsFactors = F,header=F,...)
  y = t(apply(y,1,function(x){
    x = gsub(x,pattern="^\\s+",replacement = "")
    strsplit(x,split="\\s+")[[1]]
  }))
  rownames(y) = y[,2]
  if(has_header){
    colnames(y) = y[1,]
    return(y[-1,])
  }
  return(y)
}

two_d_plot_visualize_covariate<-function(x1,x2,cov1,cov2=NULL,cuts=5,...){
  if(is.null(cov2)){cov2=cov1}
  if(is.numeric(cov1)){cov1=cut(cov1,breaks = cuts)}
  if(is.numeric(cov2)){cov1=cut(cov2,breaks = cuts)}
  cov1 = as.factor(cov1)
  cov2 = as.factor(cov2)
  cols = rainbow(length(unique(cov1)))
  names(cols) = unique(cov1)
  pchs = 1:length(unique(cov2))
  names(pchs) = unique(cov2)
  plot(x1,x2,col=cols[cov1],pch=pchs[cov2],...)
  return(list(cols,pchs))
}

cov_phe_col_to_plink_numeric_format<-function(x){
  if(is.numeric(x)){return(x)}
  num_nas = sum(is.na(x))
  v = as.numeric(as.character(x))
  if(num_nas == sum(is.na(v))){return(v)}
  x = as.numeric(as.factor(x))
  return(x)
}

from_our_sol_to_fuma_res<-function(assoc_file,bim_file,freq_file=NULL,maf = 0.005,p=1,
                                   snps_to_exclude_from_results=NULL){
  res = read.delim(assoc_file,stringsAsFactors = F)
  mafs = read.table(freq_file,stringsAsFactors = F,header=T)
  if(grepl("afreq$",freq_file)){
    colnames(mafs)[colnames(mafs)=="ID"] = "SNP"
    colnames(mafs)[colnames(mafs)=="ALT_FREQS"] = "MAF"
  }
  bim = read.delim(bim_file,stringsAsFactors = F,header=F)
  rownames(bim) = bim[,2]
  rownames(mafs) = mafs$SNP
  rownames(res) = res$ID
  selected_snps = intersect(
    rownames(res)[res$UNADJ <= p],
    rownames(mafs)[mafs$MAF >= maf]
  )
  print(paste("num selected snps in exluded set:",length(intersect(selected_snps,snps_to_exclude_from_results))))
  selected_snps = setdiff(selected_snps,snps_to_exclude_from_results)
  if(!is.element("UNADJ",set=colnames(res))){return(NULL)}
  print(table(res[selected_snps,]$UNADJ<5e-8))
  m = cbind(as.character(bim[selected_snps,1]),as.character(bim[selected_snps,4]),res[selected_snps,]$UNADJ)
  colnames(m) = c("chromosome","position","P-value")
  return(m)
}

# This function takes a vector of size 4
# The first two entries are from the first bim files and
# the second two are from the second bim file
# The function returns 1 if the snp info are the same (after sort)
# -1 if the snp is flipped
# 2 if the snp is deletion (D) or insertion (I) and ...
# 0 otherwise.
# 0 can occur, for example if we see three alleles overall
# which makes the snp problematic, especially for analysis with plink
check_bims_snp_info<-function(x){
  x1 = sort(x[1:2])
  x2 = sort(x[3:4])
  if(all(x1==x2)){return(1)}
  x1_1 = sapply(x1,rev_nucleotide)
  x1_1 = sort(x1_1)
  if(all(x2==x1_1)){return(-1)}
  return(0)
}
rev_nucleotide<-function(x){
  x = toupper(x)
  if(x=="A"){return("T")}
  if(x=="C"){return("G")}
  if(x=="G"){return("C")}
  if(x=="T"){return("A")}
  return(x)
}

# get_gwas_command_using_plink2<-function(job_dir,bfile_name){
#   curr_cmd = paste(paste(job_dir,"plink2",sep=""),
#                    "--bfile",paste(job_dir,bfile_name,sep=''),
#                    "--logistic hide-covar firth-fallback",
#                    paste("--pheno",pheno_file),
#                    paste("--pheno-name ExerciseGroup"),
#                    "--allow-no-sex",
#                    "--1",
#                    paste("--covar",covar_file),
#                    "--covar-name sex,Batch,PC1,PC2,PC3,PC4,PC5,PC6",
#                    "--adjust",
#                    "--out",paste(job_dir,"genepool_controls_simple_linear_wo_age",sep=''))
# }

create_fuma_files_for_fir<-function(dir_path,bim_file,frq_file,maf=0.005,p=0.05,
                                    snps_to_exclude_from_results=NULL){
  res_files = list.files(dir_path)
  res_files = res_files[grepl("adjusted$",res_files)]
  newdir = paste(dir_path,"fuma/",sep="")
  system(paste("mkdir",newdir))
  for (f in res_files){
    res = read.delim(paste(dir_path,f,sep=''),stringsAsFactors = F)
    res_fuma = from_our_sol_to_fuma_res(paste(dir_path,f,sep=""),bim_file,frq_file,maf = maf,p=p,
                                        snps_to_exclude_from_results=snps_to_exclude_from_results)
    write.table(res_fuma,file= paste(newdir,f,sep=""),
                row.names = F,col.names = T,quote = F,sep=" ")
  }
  
  mafs = read.table(frq_file,stringsAsFactors = F,header=T)
  mafs=mafs[,c(1,2,5,6)]
  nsample = (mafs[,3]^2 + 2*mafs[,3])*(mafs[,4]/2)
  mafs = cbind(mafs,nsample)
  write.table(mafs,file= paste(newdir,"mafs.txt",sep=""),
              row.names = F,col.names = T,quote = F,sep=" ")
  return(NULL)
}

# PCA plots
two_d_plot_visualize_covariate<-function(x1,x2,cov1,cov2=NULL,cuts=5,...){
  if(is.null(cov2)){cov2=cov1}
  n1 = length(unique(cov1))
  n2 = length(unique(cov2))
  if(is.numeric(cov1) && n1 > 50){cov1=cut(cov1,breaks = cuts)}
  if(is.numeric(cov2) && n2 > 50){cov1=cut(cov2,breaks = cuts)}
  cov1 = as.factor(cov1)
  cov2 = as.factor(cov2)
  cols = rainbow(length(unique(cov1)))
  names(cols) = as.character(unique(cov1))
  cols = cols[!is.na(names(cols))]
  pchs = 1:length(unique(cov2))
  names(pchs) = as.character(unique(cov2))
  pchs = pchs[!is.na(names(pchs))]
  plot(x1,x2,col=cols[as.character(cov1)],pch=pchs[as.character(cov2)],...)
  return(list(cols,pchs))
}
run_hclust<-function(pc_x,k,dd=NULL,h=NULL){
  if(is.null(dd)){dd = dist(pc_x,method="manhattan")}
  if(is.null(h)){h = hclust(dd,method = "complete")}
  clust = cutree(h,k=k)
  return(clust)
}
wss_score <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}
tot_wss_hluct <- function(k, hc, x) {
  cl <- cutree(hc, k)
  spl <- split(x, cl)
  wss <- sum(sapply(spl, wss_score))
  wss
}
read_pca_res<-function(path){
  pca1 = read.table(path,stringsAsFactors = F)
  r = pca1[,2]
  rownames(pca1)=r
  pca1 = pca1[,-c(1:2)]
  pca1 = as.matrix(pca1)
  if(!mode(pca1)=="numeric"){
    print("ERROR: pca matrix is not numeric");return(NULL)
  }
  colnames(pca1) = paste("PC",1:ncol(pca1),sep="")
  return(pca1)
}


################################################################################
# Sept 2018
# Create the reduced files
extract_snps_using_plink<-function(bfile,snps,out_path,snpfile,newbedfile,ref_allale_line = "",
                                  batch_script_func=get_sh_default_prefix,...){
  # create the snps file
  write.table(t(t(as.character(snps))),
              file=paste(out_path,snpfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F)
  err_path = paste(out_path,"reduce_snps",snpfile,".err",sep="")
  log_path = paste(out_path,"reduce_snps",snpfile,".log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--extract",paste(out_path,snpfile,".txt",sep=''),
                   ref_allale_line,
                   "--freq --make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"reduce_snps",snpfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
exclude_snps_using_plink<-function(bfile,snps,out_path,snpfile,newbedfile,ref_allale_line = "",
                                   batch_script_func=get_sh_default_prefix,...){
  # create the snps file
  write.table(t(t(as.character(snps))),
              file=paste(out_path,snpfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F)
  err_path = paste(out_path,"exclude_snps",snpfile,".err",sep="")
  log_path = paste(out_path,"exclude_snps",snpfile,".log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--exclude",paste(out_path,snpfile,".txt",sep=''),
                   ref_allale_line,
                   "--freq --make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"exclude_snps",snpfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
remove_subjects_using_plink<-function(bfile,subjs,out_path,subjfile,newbedfile,ref_allale_line="",
                                      batch_script_func=get_sh_default_prefix,...){
  # create the subjects file
  write.table(subjs,file=paste(out_path,subjfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F,sep="\t")
  err_path = paste(out_path,"reduce_subjs",subjfile,".err",sep="")
  log_path = paste(out_path,"reduce_subjs",subjfile,".log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--remove",paste(out_path,subjfile,".txt",sep=''),
                   ref_allale_line,
                   "--freq --make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"reduce_subjs",subjfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
keep_subjects_using_plink<-function(bfile,subjs,out_path,subjfile,newbedfile,
                                      batch_script_func=get_sh_default_prefix,...){
  # create the subjects file
  write.table(subjs,file=paste(out_path,subjfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F,sep="\t")
  err_path = paste(out_path,"reduce_subjs",subjfile,".err",sep="")
  log_path = paste(out_path,"reduce_subjs",subjfile,".log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--keep",paste(out_path,subjfile,".txt",sep=''),
                   "--freq --make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"reduce_subjs",subjfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
flip_snps_using_plink<-function(bfile,snps,out_path,snpfile,newbedfile,
                                   batch_script_func=get_sh_default_prefix,...){
  # create the snps file
  write.table(t(t(as.character(snps))),
              file=paste(out_path,snpfile,".txt",sep=''),
              row.names = F,col.names = F,quote = F)
  err_path = paste(out_path,"flip_snps_",snpfile,".err",sep="")
  log_path = paste(out_path,"flip_snps_",snpfile,".log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--flip",paste(out_path,snpfile,".txt",sep=''),
                   "--make-bed --out",paste(out_path,newbedfile,sep=''))
  curr_sh_file = paste(out_path,"flip_snps_",snpfile,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}

compute_pc_vs_binary_variable_association_p<-function(pc,y,test=wilcox.test){
  x1 = pc[y==y[1]]
  x2 = pc[y!=y[1]]
  return(test(x1,x2)$p.value)
}

compute_pc_vs_discrete_variable_association_p<-function(pc,y,cuts=5){
  if(is.null(cuts)){return(kruskal.test(pc,g=as.factor(y))$p.value)}
  pc_d = cut(pc,breaks=cuts)
  tb = table(pc_d,y)
  tb = tb[rowSums(tb)>0,]
  tb = tb[,colSums(tb)>0]
  return(chisq.test(tb)$p.value)
}

check_if_bim_is_sorted<-function(bimfile){
  d = read.table(bimfile,stringsAsFactors = F)
  d = d[d[,1]!="0",]
  ord = order(d[,1],d[,4])
  return(cor(1:length(ord),ord)>0.95)
}

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

flip_nuc<-function(x){
  if(x=="T"){return("A")}
  if(x=="A"){return("T")}
  if(x=="G"){return("C")}
  if(x=="C"){return("G")}
  if(x=="0"){return("0")}
}
flip_snp_info<-function(x){
  return(sapply(x,flip_nuc))
}
get_num_alleles<-function(x){
  return(length(setdiff(unique(x),"0")))
}

run_check_bim_analysis<-function(curr_dir,bedfile,freqfile,
                                 ref="-1000g", pop = "EUR",t_val=0.3,
                                 ref_file = "/home/users/davidama/apps/check_bim/1000GP_Phase3_combined.legend"){
  system(paste("mkdir",curr_dir))
  setwd(curr_dir)
  err_path = paste(curr_dir,"run_check_bim.err",sep="")
  log_path = paste(curr_dir,"run_check_bim.log",sep="")
  system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",curr_dir))
  # For 1000G-based analysis
  curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                   "-b", bedfile,
                   "-f", freqfile,ref,
                   "-p",pop,
                   "-t",t_val,
                   "-r",ref_file)
  curr_sh_file = "run_check_bim.sh"
  print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
                                                           mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
  system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
}
run_check_bim_output_script<-function(curr_dir,bfile_short,bfile_full){
  setwd(curr_dir)
  system(paste("less ",curr_dir,"Run-plink.sh | grep TEMP > ",curr_dir,"Run-plink2.sh",sep=""))
  run_sh_lines = readLines(paste(curr_dir,"Run-plink2.sh",sep=""))
  run_sh_lines[1] = gsub(paste("plink --bfile",bfile_short),
                         paste("plink --bfile ",bfile_full,sep=""),run_sh_lines[1])
  run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "")
  err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
  log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
  curr_sh_file = "run_check_bim_update.sh"
  print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),run_sh_lines)
  system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
}


