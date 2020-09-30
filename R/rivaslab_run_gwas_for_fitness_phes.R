library(data.table)

out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/gwas/"

rivaslab_pheno_codes_file = "/home/users/davidama/repos/ukbb-tools/05_gbe/phenotype_info.tsv"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=T)
rownames(rivaslab_codes) = rivaslab_codes[,1]

# define the self-reported activity phenotype
activity_phenotypes = c("INI884","INI894","INI904","INI914")
exercise_phenotypes = c("RestingHR" = "INI102")

ethnicity_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/sqc/population_stratification/"
ethnicity = "white_british.phe"
eth_files = list.files(ethnicity_path)
eth_file = eth_files[grepl(ethnicity,eth_files)]
eth_file = paste(ethnicity_path,eth_file,sep='')
# Get the ethnicities - subject subset to use in the current analysis
eth_subjects = fread(eth_file,header=F,stringsAsFactors = F,data.table=F)
eth_subjects = as.character(eth_subjects[,1])
print("loaded current ethnicity subjects, number of subject ids is:")
print(length(eth_subjects))

# Read the pheno
master_phe_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phewas/resources/master.phe"
phe = fread(master_phe_path,header = T,stringsAsFactors = F,data.table = F)
rownames(phe) = as.character(phe[,1])
phe = phe[intersect(eth_subjects,rownames(phe)),]
print("Dim of master phe")
print(dim(phe))
activity_phenotypes = activity_phenotypes[activity_phenotypes %in% colnames(phe)]
exercise_phenotypes = exercise_phenotypes[exercise_phenotypes %in% colnames(phe)]

# PCA of the physical activity scores
act_phe = phe[,activity_phenotypes]
act_phe = act_phe[apply(act_phe>0,1,all),]
act_pca = prcomp(act_phe,retx = T)
vars = act_pca$sdev^2
cumsum(vars/sum(vars)) # 79% of the variability in PC1 - Sept 2020
# corr between PC1 and the other scores
pc1_corrs = c(cor(act_pca$x[,1],act_phe)) # negative correlation with reported activity - Sept 2020
names(pc1_corrs) = rivaslab_codes[colnames(act_phe),2]
# make sure we are positively correlated with activity
PC1 = act_pca$x[,1]
if(sum(pc1_corrs<0)/length(pc1_corrs) >0.5){PC1 = -PC1}

# add the resting HR
restHR = phe[,"INI102"]
names(restHR) = rownames(phe)
restHR = restHR[restHR > 0]
# intersect
shared = intersect(names(restHR),names(PC1))
restHR = restHR[shared]
PC1 = PC1[shared]
# Sept 2020:
# > cor(restHR,PC1,method="spearman")
#[1] -0.033 (mild but sig)

# high values mean high activity
PC1_zs = (PC1 - mean(PC1))/sd(PC1)
# high values mean "high fitness"
restHR_zs = (restHR - mean(restHR))/sd(restHR)
restHR_zs = -1*restHR_zs
# Get the diff
activity_minus_resthr = PC1_zs - restHR_zs

# print the final phe file
curr_pcs = paste("PC",1:10,sep="")
covariate_names = c("#FID","IID","sex","age","Array",curr_pcs)
all_covs = cbind(phe[names(activity_minus_resthr),covariate_names],activity_minus_resthr)
write.table(all_covs,file=paste(out_path,"activity_minus_resthr.phe",sep=""),
                sep=" ",row.names = F,col.names = T,quote = F)

# Genetics information in a bed file
pgen_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/imp/pgen/"
all_bed_files = list.files(pgen_path)
all_bed_files = all_bed_files[grepl("bed$",all_bed_files)]  
all_bed_files = paste(pgen_path,all_bed_files,sep="")







# Ignore the rest of the file for now (Sept 2020)
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


