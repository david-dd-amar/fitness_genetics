
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/gwas/"

rivaslab_pheno_codes_file = "/home/users/davidama/repos/ukbb-tools/05_gbe/phenotype_info.tsv"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=T)
rownames(rivaslab_codes) = rivaslab_codes[,1]
activity_phenotypes = rivaslab_codes[
  grepl("activity",rivaslab_codes[,2],ignore.case = T) & 
    !grepl("anxi",rivaslab_codes[,2],ignore.case = T) &
    !grepl("depression",rivaslab_codes[,2],ignore.case = T),1]
activity_phenotypes = union(activity_phenotypes,
                            c(  "VigActivityN" = "INI1003040",
                                "ModActivityN" = "INI1003039"))
tested_phenotypes = c(
  "Exercise_HR" = "INI20005983",
  "MaxWD" = "INI30005983",
  "Recovery" = "INI40005983",
  "RestingHR" = "INI102"
)

ethnicity_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/sqc/population_stratification/"
ethnicity = "white_british.phe"
eth_files = list.files(ethnicity_path)
eth_file = eth_files[grepl(ethnicity,eth_files)]
eth_file = paste(ethnicity_path,eth_file,sep='')
# Get the ethnicities - subject subset to use in the current analysis
eth_subjects = read.delim(eth_file,header=F,stringsAsFactors = F)
eth_subjects = as.character(eth_subjects[,1])
print("loaded current ethnicity subjects, number of subject ids is:")
print(length(eth_subjects))

# Genetics information in a bed file
pgen_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/imp/pgen/"
all_bed_files = list.files(pgen_path)
all_bed_files = all_bed_files[grepl("bed$",all_bed_files)]  
all_bed_files = paste(pgen_path,all_bed_files,sep="")

# Read the pheno
master_phe_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phewas/resources/master.phe"
library(data.table,lib="~/R/packages/")
phe = fread(master_phe_path,header = T,stringsAsFactors = F,data.table = F)
tested_phenotypes[is.element(tested_phenotypes,set=colnames(phe))]
rownames(phe) = as.character(phe[,1])
phe = phe[intersect(eth_subjects,rownames(phe)),]
print("Dim of master phe")
print(dim(phe))

# Add missing phenotypes to the phe table
additional_phe_paths = "/oak/stanford/groups/mrivas/ukbb/24983/phenotypedata/extras/physical_activity/phe/"
additional_phe_files = list.files(additional_phe_paths)
additional_phe = c()
for(phe_code in c(activity_phenotypes,tested_phenotypes)){
  currfile = additional_phe_files[grepl(phe_code,additional_phe_files)]
  if(length(currfile)==0){next}
  print(paste(phe_code,"was found in the extras dir"))
  curr_v = fread(paste(additional_phe_paths,currfile,sep=""),data.table = F,header = F,stringsAsFactors = F)
  rownames(curr_v) = as.character(curr_v[,1])
  newv = rep(NA,nrow(phe))
  names(newv) = rownames(phe)
  currnames = intersect(rownames(curr_v),rownames(phe))
  newv[currnames] = curr_v[currnames,3]
  print(newv[currnames][1:10])
  additional_phe = cbind(additional_phe,newv)
  rownames(additional_phe) = names(newv)
  colnames(additional_phe)[ncol(additional_phe)] = phe_code
}
additional_phe[additional_phe==-9] = NA
phe = cbind(phe,additional_phe)

# PCA of the physical activity scores
act_phe = phe[,intersect(colnames(phe),activity_phenotypes)]
act_pca = prcomp(act_phe,retx = T)
vars = act_pca$sdev^2
cumsum(vars/sum(vars)) # 78% of the variability in PC1
# corr between PC1 and the other scores
pc1_corrs = c(cor(act_pca$x[,1],act_phe)) # positive correlation
names(pc1_corrs) = rivaslab_codes[colnames(act_phe),2]
# make sure we are positively correlated with activity
PC1 = act_pca$x[,1]
if(sum(pc1_corrs<0)/length(pc1_corrs) >0.5){PC1 = -PC1}

ExHR = additional_phe[,tested_phenotypes[1]]
rHR = phe[,tested_phenotypes[4]]
normExHR = ExHR^2 / rHR
scoresMat = cbind(ExHR,rHR,normExHR,PC1)
scoresMat = scoresMat[!apply(is.na(scoresMat),1,any),]
cor(scoresMat)

tb = table(normExHR < quantile(normExHR,na.rm = T,probs = 0.1),
      PC1 < quantile(PC1,na.rm = T,probs = 0.1))
# tb = table(ExHR > quantile(ExHR,na.rm = T,probs = 0.9),
#            rHR > quantile(rHR,na.rm = T,probs = 0.9))
fisher.test(tb)

# Create folders with a phe file and a covariate file
PCs = c(5,10)
for(npc in PCs){
  curr_pcs = paste("PC",1:npc,sep="")
  curr_covariate_names = c("sex","age","Array",curr_pcs)
  for(pheno_name in names(tested_phenotypes)){
    pheno_code = tested_phenotypes[pheno_name]
    if(!is.element(pheno_code,set=colnames(phe))){next}
    curr_path = paste(out_path,pheno_name,"_PCs",npc,"/",sep="")
    try(system(paste("mkdir",curr_path)))
    # create the phe file and the covars file
    y = phe[,pheno_code]
    x = phe[,curr_covariate_names]
    x = as.matrix(x)
    rownames(x) = rownames(phe)
    names(y) = rownames(phe)
    # If the phenotype is numeric, regress out the covariats
    newy = y
    covs = NULL
    if(mode(y)=="numeric" && length(unique(y))>2){
      lm = lm(y~x)
      residuals = lm$residuals
      newy = residuals
    }
    else{
      covs = cbind(rownames(phe),rownames(phe),x)
      colnames(covs) = c("#FID","IID",curr_covariate_names)
    }
    newy = cbind(names(newy),names(newy),newy)
    colnames(newy) = c("#FID","IID",pheno_name)
    write.table(newy,file=paste(curr_path,pheno_name,".phe",sep=""),
                sep=" ",row.names = F,col.names = T,quote = F)
    if(! is.null(covs)){
      write.table(covs,file=paste(curr_path,"covs.txt",sep=""),
                  sep=" ",row.names = F,col.names = T,quote = F)
    }
  }
}

###############################################################################
###############################################################################
###############################################################################
# Useful functions for our analyses below
print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}

get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",
              Ncpu,mem_size,plink_pkg = "plink/2.0a1"){
  partition_line = "#SBATCH --partition=euan,mrivas,normal,owners"
  if(mem_size>128000){
    partition_line = "#SBATCH --partition=bigmem,euan,mrivas"
  }
  return(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH --time=18:00:00",
      partition_line,
      "#SBATCH --nodes=1",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "#SBATCH -x sh-113-15",
      "",
      "module load biology",
      paste("module load",plink_pkg)
    )
  )
}
###############################################################################
###############################################################################
###############################################################################

# Run the GWAS  
PCs = c(5,10)
for(npc in PCs){
  curr_pcs = paste("PC",1:npc,sep="")
  curr_covariate_names = c("sex","age","Array",curr_pcs)
  for(pheno_name in names(tested_phenotypes)){
    pheno_code = tested_phenotypes[pheno_name]
    if(!is.element(pheno_code,set=colnames(phe))){next}
    curr_path = paste(out_path,pheno_name,"_PCs",npc,"/",sep="")
    
    regr_cmd_line = "--linear hide-covar"
    covars_line1 = ""; covars_line2 = ""
    if(is.element("covs.txt",set=list.files("curr_path"))){
      regr_cmd_line = "--logistic hide-covar firth-fallback"
      covars_line1 = paste("--covar",paste(curr_path,"covs.txt",sep=""))
      covars_line2 = paste("--covar-name",paste(curr_covariate_names,collapse=","))
    }
    
    for(chr in all_bed_files){
      bed = gsub(".bed$","",chr)
      chr_regx = regexpr("chr\\d+",chr)
      taskname = substr(chr,chr_regx,chr_regx + attr(chr_regx,"match.length")-1)
      err_path = paste(curr_path,taskname,"_tmp.err",sep="")
      log_path = paste(curr_path,taskname,"_tmp.log",sep="")
      curr_cmd = paste("plink2",
                       "--bfile",bed,
                       regr_cmd_line,
                       "--pheno", paste(curr_path,pheno_name,".phe",sep=""),
                       "--pheno-name", pheno_name, 
                       covars_line1,covars_line2,
                       "--threads 4",
                       "--maf 0.01",
                       "--adjust",
                       "--out",paste(curr_path,pheno_name,"_",taskname,sep=''))
      curr_sh_file = paste(curr_path,taskname,"_tmp.sh",sep="")
      print_sh_file(curr_sh_file,
                    get_sh_prefix_one_node_specify_cpu_and_mem(
                      err_path,log_path,Ncpu=4,mem_size=32000),
                    curr_cmd)
      system(paste("sbatch",curr_sh_file))
    }
  }
}


