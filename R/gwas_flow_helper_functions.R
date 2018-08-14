####################################################################################################
####################################################################################################
####################################################################################################
######## Some useful functions
get_sh_default_prefix<-function(err="",log=""){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=12:00:00",
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=2",
      "#SBATCH --mem=32000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}

# plink2: plink/2.0a1
get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",plink_pkg = "plink/1.90b5.3",Ncpu,mem_size){
  partition_line = "#SBATCH --partition=euan,mrivas,normal,owners"
  if(mem_size>32000){
    partition_line = "#SBATCH --partition=bigmem,euan,mrivas"
  }
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=6:00:00",
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

get_sh_prefix_bigmem<-function(err="",log="",plink_pkg = "plink/1.90b5.3",Ncpu=1,mem_size){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=6:00:00",
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
  #new_jobs = rownames(jobs)[jobs[,5]=="RUNNING" | jobs[,5]=="PENDING"]
  return(jobs)
}
get_job_id<-function(x){return(x[1])}
wait_for_job<-function(jobs_before=NULL,waittime=6,max_wait=6000){
  Sys.sleep(waittime)
  curr_jobs = get_my_jobs()
  new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
  if(length(new_jobs)==0){return(NULL)}
  print(paste("new added jobs are: ",new_jobs))
  i = 1
  while(T){
    Sys.sleep(waittime)
    curr_jobs = get_my_jobs()
    new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
    if(length(new_jobs)==0){return(NULL)}
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

from_our_sol_to_fuma_res<-function(assoc_file,bim_file,freq_file=NULL,maf = 0.001,p=1){
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

create_fuma_files_for_fir<-function(dir_path,bim_file,frq_file,maf=0.001){
  res_files = list.files(dir_path)
  res_files = res_files[grepl("adjusted$",res_files)]
  newdir = paste(dir_path,"fuma/",sep="")
  system(paste("mkdir",newdir))
  for (f in res_files){
    res = read.delim(paste(dir_path,f,sep=''),stringsAsFactors = F)
    res_fuma = from_our_sol_to_fuma_res(paste(dir_path,f,sep=""),bim_file,frq_file,maf = maf)
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
  if(is.numeric(cov1) && n1 > 10){cov1=cut(cov1,breaks = cuts)}
  if(is.numeric(cov2) && n2 > 10){cov1=cut(cov2,breaks = cuts)}
  cov1 = as.factor(cov1)
  cov2 = as.factor(cov2)
  cols = rainbow(length(unique(cov1)))
  names(cols) = unique(cov1)
  cols = cols[!is.na(names(cols))]
  pchs = 1:length(unique(cov2))
  names(pchs) = unique(cov2)
  pchs = pchs[!is.na(names(pchs))]
  plot(x1,x2,col=cols[cov1],pch=pchs[cov2],...)
  return(list(cols,pchs))
}
