# We do a quick QC analysis of selected chromosomes
# Our goal is to illustrate that intra-eur population structure
# is crucial for the analysis of our dataset
# The analysis here is straightforward. We select data from three chromosomes,
# take only snps that are not palindromic and appear in both datasets.
# We then merge the datasets and do PCA.
# All analyses are done using plink.

# ASSUMPTIONS: dataset1 is 1000g, dataset2 is mega-related dataset AFTER running
# check_bim. Thus, most alleles should be matched by now.

# Load auxiliary functions
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)
# Define the datasets to be analyzed
dataset1 = "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/"
# dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_eu_imp_wo_jhu/with_ukbb/"
# dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/by_chr/"
dataset2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/by_chr/"
chrs = paste("chr",c(1:22),sep="")
# define the output path
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/qc/1000g_pop/"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/qc/1000g_vs_direct_mega/"
system(paste("mkdir",out_path))

# mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/integrated_sample_metadata_and_covariates.phe"
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"


is_snp_paly<-function(x){
	return(x=="AT" || x=="TA" || x=="CG" || x=="GC")
}

# Analyze the datasets and prepare for a merge
for(chr in chrs){
	bim1 = read.table(paste(dataset1,chr,".bim",sep=""),stringsAsFactors=F)
	bim2 = read.table(paste(dataset2,chr,".bim",sep=""),stringsAsFactors=F)
	raw_id_intersect = intersect(bim1[,2],bim2[,2])
	bim1_dups = names(which(table(bim1[,2])>1))
	bim2_dups = names(which(table(bim2[,2])>1))
	raw_id_intersect = setdiff(raw_id_intersect,bim1_dups)
	raw_id_intersect = setdiff(raw_id_intersect,bim2_dups)
	bim1 = bim1[is.element(bim1[,2],set=raw_id_intersect),]
	bim2 = bim2[is.element(bim2[,2],set=raw_id_intersect),]
	rownames(bim1) = bim1[,2]
	rownames(bim2) = bim2[,2]
	bim1 = bim1[raw_id_intersect,]
	bim2 = bim2[raw_id_intersect,]
	all(bim1[,2]==bim2[,2]) # sanity check
	
	# Exclude palindromic SNPs
	alleles2 = paste(bim2[,5],bim2[,6],sep="")
	pali_snps = raw_id_intersect[sapply(alleles2,is_snp_paly)]
	raw_id_intersect = setdiff(raw_id_intersect,pali_snps)
	bim1 = bim1[raw_id_intersect,]
	bim2 = bim2[raw_id_intersect,]
	
	alleles1 = paste(bim1[,5],bim1[,6],sep="")
	rev_alleles1 = paste(bim1[,6],bim1[,5],sep="")
	alleles2 = paste(bim2[,5],bim2[,6],sep="")
	table(rev_alleles1 == alleles2) # almost all should be TRUE
	
	# NOTE: if we want to add the ability to handle snp flip for the merge
	# then we need another file with a list of snps to flip in one of the files
	# this will also include matching the alleles to make sure we get the correct
	# ones for the flip
	snps_to_take = raw_id_intersect[rev_alleles1 == alleles2 | alleles1 == alleles2]
	curr_snps_file = paste(out_path,chr,"_curr_snps.txt",sep="")
	write.table(t(t(snps_to_take)),file=curr_snps_file,row.names=F,col.names=F,quote=F)
	print(paste("chromosome:",chr,"number of snps:",length(snps_to_take)))
	# create a bed from dataset1
	curr_name = paste("dataset1_",chr,"_make_bed",sep='')
  	curr_cmd = paste("plink --bfile",paste(dataset1,chr,sep=""),
                   "--extract",curr_snps_file,
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"dataset1_",chr,sep='')
  	)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
    
    # create a bed from dataset2
    curr_name = paste("dataset2_",chr,"_make_bed",sep='')
  	curr_cmd = paste("plink --bfile",paste(dataset2,chr,sep=""),
                   "--extract",curr_snps_file,
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"dataset2_",chr,sep='')
  	)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(240)

# Compute freqs
for(chr in chrs){
	f1 = paste(out_path,"dataset1_",chr,sep='')
	f2 = paste(out_path,"dataset2_",chr,sep='')
	
	curr_name = paste("dataset1_",chr,"_frq",sep='')
  	curr_cmd = paste("plink --bfile",f1,
                   "--freq","--threads 4",
                   "--out",f1)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
    
	curr_name = paste("dataset2_",chr,"_frq",sep='')
  	curr_cmd = paste("plink --bfile",f2,
                   "--freq","--threads 4",
                   "--out",f2)
  	run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)	
}

# Do we have SNPs that are essentially zero maf in one dataset and 
# high in the other?
maf_thr1 = 0.001;maf_thr2 = 0.05
problematic_snps = c()
for(chr in chrs){
  f1 = paste(out_path,"dataset1_",chr,".frq",sep='')
  f2 = paste(out_path,"dataset2_",chr,".frq",sep='')
  freqs1 = read.table(f1,header=T)
  freqs2 = read.table(f2,header=T)
  rownames(freqs1) = freqs1$SNP
  rownames(freqs2) = freqs2$SNP
  inds = intersect(freqs1$SNP,freqs2$SNP)
  x1 = freqs1[inds,"MAF"]
  x2 = freqs2[inds,"MAF"]
  problematic_snps = c(problematic_snps,
                       inds[(x1 >=maf_thr2 & x2 < maf_thr1) | (x1 < maf_thr1 & x2 >= maf_thr2)]
  )
}
length(problematic_snps)

# Merge each chromosome data
for(chr in chrs){
	f1 = paste(out_path,"dataset1_",chr,sep='')
	f2 = paste(out_path,"dataset2_",chr,sep='')
	curr_name = paste("merge_",chr,sep='')
	curr_cmd = paste("plink --bfile",f1,
                   "--bmerge",f2,"--threads 8",
                   "--make-bed --out",paste(out_path,"merged_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Run LD prunning on each new file
for(chr in chrs){
	curr_name = paste("ld_",chr,sep='')
	curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--threads 8",
                   "--maf 0.05",
                   "--indep-pairwise 500 10",0.1,
                   "--out",paste(out_path,"merged_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Compare pruned results to the problematic SNPs above
all_selected_snps_for_pca = c()
for(chr in chrs){
  curr_snps = read.table(paste(out_path,"merged_",chr,".prune.in",sep=''),stringsAsFactors = F)[,1]
  all_selected_snps_for_pca = c(all_selected_snps_for_pca,curr_snps)
}
print(length(all_selected_snps_for_pca))
print(intersect(all_selected_snps_for_pca,problematic_snps))

# Extract LD-pruned beds
for(chr in chrs){
	curr_name = paste("ld_extract_",chr,sep='')
	curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--threads 8",
                   "--extract", paste(out_path,"merged_",chr,".prune.in",sep=''),
                   "--make-bed --out",paste(out_path,"merged_ld_pruned_",chr,sep='')
	)
	run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)
}

# Merge all resulting bed files into a single one
all_bfiles = paste(out_path,"merged_ld_pruned_",chrs,sep='')
write.table(t(t(all_bfiles[-1])),file=paste(out_path,"all_bfiles.txt",sep=""),
	row.names=F,col.names=F,quote=F)
curr_name = paste("merge_beds",sep='')
curr_cmd = paste("plink --bfile",all_bfiles[1],
                   "--merge-list",paste(out_path,"all_bfiles.txt",sep=""),
                   "--threads 4",
                   "--make-bed --out",paste(out_path,"merged_dataset",sep='')
)
run_plink_command(curr_cmd,out_path,curr_name,
                get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=64000)

# Run PCA 
curr_name = "run_pca"
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
        	    "--threads 16",
                "--pca 40",
                "--out",paste(out_path,curr_name,sep='')
)
run_plink_command(curr_cmd,out_path,curr_name,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)

#####################################################################################
#####################################################################################
#####################################################################################
pcax = read_pca_res(paste(out_path,"run_pca.eigenvec",sep=""))
d2 = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/1000g/integrated_call_samples_v3.20130502.ALL.panel"
                ,stringsAsFactors=F,header=T,row.names = 1)
mega_covars = read.table(mega_covars_path,header=T)
rownames(mega_covars) = mega_covars[,2]

cohorts2 = c(rownames(mega_covars),d2[,1])
names(cohorts2) = cohorts2
cohorts2[rownames(mega_covars)] = mega_covars$Cohort
cohorts2[rownames(d2)] = d2[,1]


# Define EUR-descendants manually
yy = pcax[,1] < 0 & pcax[,2]< 0
table(yy)
table(yy,cohorts2[names(yy)])

# For Cooper and Elite, define the subjects are rerun the PCA
our_subjects = names(yy)[yy]
curr_fam = read.table(paste(out_path,"merged_dataset.fam",sep=''))
curr_fam = curr_fam[is.element(curr_fam[,2],set=our_subjects),]
dim(curr_fam)

curr_path = paste(out_path,"gwas_eu/",sep="")
system(paste("mkdir",curr_path))
write.table(curr_fam,file=paste(curr_path,"our_subjects.txt",sep=""),
            row.names = F,col.names = F,quote=F,sep=" ")

# Rerun the PCA, add curr maf filter
# Run PCA 
curr_name = "rerun_pca"
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
                 "--threads 16",
                 "--keep",paste(curr_path,"our_subjects.txt",sep=""),
                 "--maf 0.05",
                 "--pca 40",
                 "--out",paste(curr_path,curr_name,sep='')
)
run_plink_command(curr_cmd,curr_path,curr_name,
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=16,mem_size=64000)
wait_for_job()

# put the covars in a single file
new_pcax  = read_pca_res(paste(curr_path,"rerun_pca.eigenvec",sep=""))
rownames(curr_fam) = curr_fam[,2]
sex = rep(NA,nrow(new_pcax))
names(sex) = rownames(new_pcax)
inds1 = intersect(names(sex),rownames(mega_covars))
sex[inds1] = mega_covars[inds1,"sex"]
inds2 = intersect(rownames(d2),names(sex))
sex[inds2] = d2[inds2,3]
table(sex)
sex[sex=="female"]="2"
sex[sex=="male"]="1"
table(sex)

curr_cohorts = cohorts2[rownames(new_pcax)]
cooper_col = rep(NA,nrow(new_pcax))
names(cooper_col) = rownames(new_pcax)
cooper_col[curr_cohorts=="1"]  = 1
cooper_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
elite_col = rep(NA,nrow(new_pcax))
names(elite_col) = rownames(new_pcax)
elite_col[curr_cohorts=="2"]  = 1
elite_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2
gp_col = rep(NA,nrow(new_pcax))
names(gp_col) = rownames(gp_col)
gp_col[curr_cohorts=="3"]  = 1
gp_col[curr_cohorts !="1" & curr_cohorts!= "2" & curr_cohorts!= "3"]  = 2

covs = cbind(curr_fam[rownames(new_pcax),],sex,new_pcax,cooper_col,elite_col,gp_col)
colnames(covs)[1:2] = c("FID","IID")
covs_file = "our_subjects_covs.phe"
write.table(covs,file=paste(curr_path,covs_file,sep=""),
            row.names = F,col.names = T,quote=F,sep=" ")
covs_file = paste(curr_path,"our_subjects_covs.phe",sep="")

# Run the GWAS
for(chr in chrs){
  curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--logistic hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name cooper_col",
                   "--covar",covs_file,
                   "--maf 0.01",
                   "--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                   "--allow-no-sex --adjust",
                   "--threads",4,
                   "--out",paste(curr_path,"cooper_gwas_res_",chr,sep="")
  )
  run_plink_command(curr_cmd,curr_path,paste("cooper_gwas_res_",chr,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
}
for(chr in chrs){
  curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--logistic hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name elite_col",
                   "--covar",covs_file,
                   "--maf 0.01",
                   "--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                   "--allow-no-sex --adjust",
                   "--threads",4,
                   "--out",paste(curr_path,"elite_gwas_res_",chr,sep="")
  )
  run_plink_command(curr_cmd,curr_path,paste("elite_gwas_res_",chr,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
}
for(chr in chrs){
  curr_cmd = paste("plink --bfile",paste(out_path,"merged_",chr,sep=''),
                   "--logistic hide-covar",
                   "--pheno",covs_file,
                   "--pheno-name gp_col",
                   "--covar",covs_file,
                   "--maf 0.01",
                   "--covar-name sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
                   "--allow-no-sex --adjust",
                   "--threads",4,
                   "--out",paste(curr_path,"gp_gwas_res_",chr,sep="")
  )
  run_plink_command(curr_cmd,curr_path,paste("gp_gwas_res_",chr,sep=""),
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
}

# Merge the results into single files
res_files = paste(curr_path,"elite_gwas_res_",chrs,".assoc.logistic",sep="")
res_file = paste(curr_path,"elite_gwas_res_all.assoc",sep="")
for(j in 1:length(res_files)){
  if(j==1){
    system(paste("less",res_files[j],">",res_file))
  }
  if(j>1){
    system(paste("less",res_files[j],"| grep -v SNP >>",res_file))
  }
}
# Merge the results into single files
res_files = paste(curr_path,"cooper_gwas_res_",chrs,".assoc.logistic",sep="")
res_file = paste(curr_path,"cooper_gwas_res_all.assoc",sep="")
for(j in 1:length(res_files)){
  if(j==1){
    system(paste("less",res_files[j],">",res_file))
  }
  if(j>1){
    system(paste("less",res_files[j],"| grep -v SNP >>",res_file))
  }
}
# GP: Merge the results into single files
res_files = paste(curr_path,"gp_gwas_res_",chrs,".assoc.logistic",sep="")
res_file = paste(curr_path,"gp_gwas_res_all.assoc",sep="")
for(j in 1:length(res_files)){
  if(j==1){
    system(paste("less",res_files[j],">",res_file))
  }
  if(j>1){
    system(paste("less",res_files[j],"| grep -v SNP >>",res_file))
  }
}
# write into fuma files
res_file = paste(curr_path,"elite_gwas_res_all.assoc",sep="")
res_file2 = paste(curr_path,"fuma_elite_gwas_res_all.assoc",sep="")
res = read.table(res_file,header=T,stringsAsFactors = F)
res_elite=res
res = res[,c("SNP","P")]
colnames(res) = c("rsID","P-value")
write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")

res_file = paste(curr_path,"cooper_gwas_res_all.assoc",sep="")
res_file2 = paste(curr_path,"fuma_cooper_gwas_res_all.assoc",sep="")
res = read.table(res_file,header=T,stringsAsFactors = F)
res_cooper=res
res = res[,c("SNP","P")]
colnames(res) = c("rsID","P-value")
write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")

res_file = paste(curr_path,"gp_gwas_res_all.assoc",sep="")
res_file2 = paste(curr_path,"fuma_gp_gwas_res_all.assoc",sep="")
res = read.table(res_file,header=T,stringsAsFactors = F)
res_gp=res
res = res[,c("SNP","P")]
colnames(res) = c("rsID","P-value")
write.table(res,file=res_file2,col.names = T,row.names = F,quote = F,sep=" ")

rownames(res_cooper) = res_cooper$SNP
rownames(res_elite) = res_elite$SNP
rownames(res_gp) = res_gp$SNP
snps = intersect(res_gp$SNP,res_elite$SNP)
x1 = res_cooper[snps,"P"]
x2 = res_elite[snps,"P"]
x2 = res_gp[snps,"P"]
table(x1<1e-5,x2<1e-5)

d = read.table(paste(curr_path,"our_subjects_covs.phe",sep=""),header=T)
pc_rocs = c()
for(j in 1:40){
  curr_inds = !is.na(d[,"elite_col"]) 
  p1 = compute_pc_vs_binary_variable_association_roc(
    pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"elite_col"]
  )
  curr_inds = !is.na(d[,"cooper_col"]) 
  p2 = compute_pc_vs_binary_variable_association_roc(
    pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"cooper_col"]
  )
  pc_rocs = rbind(pc_ps,c(p1,p2))
}
pc_qs = apply(pc_ps,2,p.adjust)
pc_inds = pc_qs < 0.01

#####################################################################################
#####################################################################################
#####################################################################################

# Take a specific PC after manual inspection and run a GWAS against it
PC = "PC6"
curr_path = paste(out_path,"gwas_",PC,"/",sep="")
system(paste("mkdir",curr_path))
pca_data = read.table(paste(out_path,"run_pca.eigenvec",sep=''))
colnames(pca_data) = c("FID","IID",paste("PC",1:40,sep=""))
write.table(pca_data,paste(curr_path,"pcs.phe",sep=""),
            row.names = F,col.names = T,quote=F,sep="\t")
curr_cmd = paste("plink --bfile",paste(out_path,"merged_dataset",sep=''),
  "--linear hide-covar",
  paste("--pheno",paste(curr_path,"pcs.phe",sep="")),
  "--pheno-name",PC,
  "--threads",8,
  "--out",paste(curr_path,"gwas_res",sep="")
)
run_plink_command(curr_cmd,curr_path,"pc_gwas_run",
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)
wait_for_job(120)
# read the results, look at the selected SNPs
res = read.table(paste(curr_path,"gwas_res.assoc.linear",sep=""),
                 stringsAsFactors = F,row.names=2,header=T)
curr_PC_snps = rownames(res)[res$P < 1e-8]
for(chr in chrs){
  f1 = paste(out_path,"dataset1_",chr,".frq",sep='')
  f2 = paste(out_path,"dataset2_",chr,".frq",sep='')
  freqs1 = read.table(f1,header=T)
  freqs2 = read.table(f2,header=T)
  rownames(freqs1) = freqs1$SNP
  rownames(freqs2) = freqs2$SNP
  inds = intersect(freqs1$SNP,freqs2$SNP)
  inds = intersect(inds,curr_PC_snps)
  break
}

# are these directly genotypes?
direct_g = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data.bim",
                      row.names = 2,stringsAsFactors = F)
intersect(curr_PC_snps,rownames(direct_g))


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Locally
setwd("/Users/David/Desktop/elite/november2018_analysis/qc_1000g_mega_direct/")
source("~/Desktop/repos/fitness_genetics/R/gwas_flow_helper_functions.R")

pcax = read_pca_res("run_pca.eigenvec")
d = read.table("/Users/David/Desktop/elite/november2018_analysis/qc_1000g/all_cohorts.phe",header=T,stringsAsFactors = F)
d = read.table("integrated_sample_metadata_and_covariates.phe",header=T,stringsAsFactors = F)
rownames(d) = d$IID
d2 = read.table("/Users/David/Desktop/elite/november2018_analysis/qc_1000g/integrated_call_samples_v3.20130502.ALL.panel"
                ,stringsAsFactors=F,header=T)
inds = intersect(rownames(d),rownames(pcax))
d = d[inds,]
dim(d)

CohortName = d$Cohort
table(CohortName)
CohortName[CohortName=="1"]="Cooper"
CohortName[CohortName=="2"]="ELITE"
d = cbind(d,CohortName)
d$CohortName = as.character(d$CohortName)

cohorts1 = c(d$CohortName,d2[,2])
names(cohorts1) = c(rownames(d),d2[,1])
cohorts2 = c(d$CohortName,d2[,3])
names(cohorts2) = c(rownames(d),d2[,1])
cohorts3 = c(d$CohortName,rep("1000g",nrow(d2)))
names(cohorts3) = c(rownames(d),d2[,1])
pcax = pcax[names(cohorts1),]

# PCA plots
inds = rownames(pcax)
# inds = names(cohorts2)[cohorts2=="elite" | cohorts2=="cooper" | cohorts2=="EUR"]
res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],cohorts2[inds],cohorts2[inds],
      main = "All filters",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,3],pcax[inds,4],cohorts2[inds],cohorts2[inds],
      main = "All filters",xlab="PC3",ylab="PC4",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,5],pcax[inds,6],cohorts2[inds],cohorts2[inds],
      main = "All filters",xlab="PC5",ylab="PC6",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,7],pcax[inds,8],cohorts2[inds],cohorts2[inds],
    main = "All filters",xlab="PC7",ylab="PC8",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,9],pcax[inds,10],cohorts2[inds],cohorts2[inds],
    main = "All filters",xlab="PC9",ylab="PC10",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,11],pcax[inds,12],cohorts2[inds],cohorts2[inds],
    main = "All filters",xlab="PC11",ylab="PC12",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

res = two_d_plot_visualize_covariate(pcax[inds,13],pcax[inds,14],cohorts2[inds],cohorts2[inds],
    main = "All filters",xlab="PC13",ylab="PC14",lwd=2)
legend(x="bottomleft",names(res[[1]]),fill = res[[1]])

# Define EUR-descendants manually
yy = pcax[,1] < 0 & pcax[,2]< -0.002
table(yy,cohorts2[names(yy)])
table(yy,cohorts1[names(yy)])
table(yy)

res = two_d_plot_visualize_covariate(pcax[inds,1],pcax[inds,2],yy[inds],yy[inds],
  main = "All filters",xlab="PC1",ylab="PC2",lwd=2)
legend(x="topright",names(res[[1]]),fill = res[[1]])


# Do the pcs perfectly predict the group?
y = as.factor(cohorts3=="ukbb")
y = as.factor(cohorts3=="elite")
y = as.factor(cohorts3=="cooper")
y = as.factor(cohorts3=="1000g")
x = pcax[,1:40]
logistic_d = data.frame(y,x)
cv_res = c()
folds = sample(rep(1:10,length(y)/10))[1:length(y)]
for(i in 1:10){
  tr = logistic_d[folds!=i,]
  te = x[folds==i,]
  tey = y[folds==i]
  m = glm(y~.,family = "binomial",data=tr)
  preds = predict(m,newdata=data.frame(logistic_d[folds==i,]),type="response")
  cv_res = rbind(cv_res,cbind(preds,tey))
}
table(cv_res[,1]>0.5,cv_res[,2])

res_cooper = read.table("fuma_cooper_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)
res_elite = read.table("fuma_elite_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)
res_gp = read.table("fuma_gp_gwas_res_all.assoc",stringsAsFactors = F,header=T,row.names = 1)

inds = intersect(rownames(res_cooper),rownames(res_elite))
qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_elite$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)

inds = intersect(inds,rownames(res_gp))
qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_elite$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)
qqplot(y=-log(res_gp$P.value,10),x=-log(runif(100000),10),pch=20,
       ylab="Sample quantiles",xlab="Theoretic quantiles");abline(0,1)

plot(y=-log(res_gp[inds,"P.value"],10),x=-log(res_elite[inds,"P.value"],10),pch=20,
       ylab="GP p-value",xlab="ELITE p-value");abline(0,1)

table(res_gp[inds,"P.value"] < 10^-5,res_elite[inds,"P.value"]<10^-5 )
na.omit(inds[res_gp[inds,"P.value"] < 10^-5])

# qqplot(y=-log(res_cooper$P.value,10),x=-log(runif(100000),10),pch=20,
#        ylab="Sample quantiles",xlab="Theoretic quantiles",
#        xlim=c(0,20),ylim=c(0,20));abline(0,1)

# library(lattice);
# x=res$P
# qqmath(~-log10(na.omit(x)),
#        distribution=function(x){-log10(qunif(1-x))}
# )
# abline(0,1,add=T)









































