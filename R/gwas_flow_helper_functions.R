####################################################################################################
####################################################################################################
####################################################################################################
######## Some useful functions
get_sh_default_prefix<-function(err="",log=""){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=24:00:00",
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=2",
      "#SBATCH --mem=64000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}

get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",Ncpu,mem_size){
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=24:00:00",
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      paste("#SBATCH --c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}

print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}

get_my_jobs<-function(){
  tmp = paste("tmp",abs(rnorm(1)),sep='')
  system(paste("squeue | grep davidama >",tmp),wait = T)
  jobs = readLines(tmp)
  system(paste("rm",tmp))
  return(jobs)
}

get_job_id<-function(x){
  x = strsplit(x,split="\\s+",perl=T)[[1]][-1]
  return(x[1])
}

wait_for_job<-function(jobs_before,waittime=30){
  Sys.sleep(waittime)
  jobs_before = sapply(jobs_before,get_job_id)
  curr_jobs = get_my_jobs()
  curr_jobs = sapply(curr_jobs,get_job_id)
  new_job = setdiff(curr_jobs,jobs_before)
  if(length(new_job)==0){return(NULL)}
  print(paste("new added job is: ",new_job))
  while(is.element(new_job,set=curr_jobs)){
    Sys.sleep(waittime)
    curr_jobs = get_my_jobs()
    curr_jobs = sapply(curr_jobs,get_job_id)
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
  x = as.numeric(as.factor(x))
  x[is.na(x)] = -9
  return(x)
}



