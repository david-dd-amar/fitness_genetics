library(data.table)

out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/gwas/"
computed_fitness_scores = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_exercise_scores.phe"

rivaslab_pheno_codes_file = "/home/users/davidama/repos/ukbb-tools/05_gbe/phenotype_info.tsv"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=T)
rownames(rivaslab_codes) = rivaslab_codes[,1]

comp_exercise_data = fread(computed_fitness_scores,stringsAsFactors = F,header=T,data.table=F)

# define the self-reported activity phenotype
activity_phenotypes = c("NdaysModerate" = "INI884","DurModerate"="INI894","NdaysVig" = "INI904", "DurVig" = "INI914")
exercise_phenotypes = c("RestingHR" = "INI102")
background_variables = c("Height" = "INI50", "Weight" = "INI21002","WaistCirc"="INI48","BMI" = "INI21001","age"="age","sex"="sex")
for(j in 1:20){
  currpc = paste0("PC",j)
  background_variables[currpc]=currpc
}
computed_exercise_scores = colnames(comp_exercise_data)[-c(1,2)]

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

# merge all datasets
phe_subset = phe[,c("#FID","IID",background_variables,activity_phenotypes,exercise_phenotypes)]
colnames(phe_subset) = c("#FID","IID",names(background_variables),names(activity_phenotypes),names(exercise_phenotypes))
comp_exercise_data_eu = comp_exercise_data[comp_exercise_data[,1] %in% phe_subset[,1],]
merged_phe = merge(phe_subset,comp_exercise_data_eu,by="IID",all=T)
all(merged_phe[,2] == merged_phe$FID,na.rm=T)
merged_phe = merged_phe[,colnames(merged_phe)!="FID"]
merged_phe = merged_phe[,c(2:1,3:ncol(merged_phe))]

write.table(merged_phe,file=paste(out_path,"all_variables_for_analysis.phe",sep=""),
                sep=" ",row.names = F,col.names = T,quote = F)


# Older analysis: attempt to look at activity and rhr only
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

# same as before but using normal z-scores
# rank - high value == high score
PC1_ranks = rank(PC1)/length(PC1)
restHR_ranks = rank(restHR)/length(PC1)
# fix 1 and 0 issues
PC1_ranks[PC1_ranks==1] = max(PC1_ranks[PC1_ranks < 1])
PC1_ranks[PC1_ranks==0] = min(PC1_ranks[PC1_ranks > 0])
restHR_ranks[restHR_ranks==1] = max(restHR_ranks[restHR_ranks < 1])
restHR_ranks[restHR_ranks==0] = min(restHR_ranks[restHR_ranks > 0])
PC1_normz = qnorm(PC1_ranks)
restHR_normz = qnorm(restHR_ranks)
# these z-scores should be positively correlated with the actual scores:
# sept 2020
#> cor(PC1_normz,PC1)
#[1] 0.8357246
#cor(restHR,restHR_normz)
#[1] 0.9909848
activity_minus_resthr_znorm = PC1_normz - restHR_normz

# print the final phe file
curr_pcs = paste("PC",1:10,sep="")
covariate_names = c("#FID","IID","sex","age","Array",curr_pcs)
all_covs = cbind(phe[names(activity_minus_resthr),covariate_names],activity_minus_resthr,activity_minus_resthr_znorm)
activityPC1 = PC1
all_covs = cbind(all_covs,restHR,activityPC1)
write.table(all_covs,file=paste(out_path,"activity_minus_resthr.phe",sep=""),
                sep=" ",row.names = F,col.names = T,quote = F)

# Analyze the exercise test data
###############################################################################
###############################################################################
###############################################################################

# Compare GWAS results
###############################################################################
###############################################################################
###############################################################################
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/gwas/"
setwd(out_path)
library(data.table)

gfiles = c(
  "exercise.activity_minus_resthr_znorm.glm.linear",
  "rhr.restHR.glm.linear",
  "activityPC1.activityPC1.glm.linear",
  "exercise.activity_minus_resthr.glm.linear"
)

gdata = list()
for(gfile in gfiles){
   gdata[[gfile]] = fread(gfile,data.table=F,header=T,stringsAsFactors=F)
}
for(gfile in gfiles){
   rownames(gdata[[gfile]]) = gdata[[gfile]][,"ID"]
}

sapply(gdata,names)
pthr = 1e-08
variant_sets = lapply(gdata,function(x,pthr)x[x$P < pthr,"ID"],pthr=pthr)
sapply(variant_sets,length)
length(setdiff(variant_sets[[1]],variant_sets[[2]]))
length(setdiff(variant_sets[[4]],variant_sets[[2]]))
set1 = setdiff(variant_sets[[1]],variant_sets[[2]])
gdata[[2]][set1,"P"]

cor(gdata[[1]]$P,gdata[[2]]$P)
cor(gdata[[1]]$BETA,gdata[[2]]$BETA)








# Ignore the rest of the file for now (Sept 2020)
###############################################################################
###############################################################################
###############################################################################

# Genetics information in a bed file
pgen_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/imp/pgen/"
all_bed_files = list.files(pgen_path)
all_bed_files = all_bed_files[grepl("bed$",all_bed_files)]  
all_bed_files = paste(pgen_path,all_bed_files,sep="")


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


