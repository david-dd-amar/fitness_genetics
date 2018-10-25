
## Look at rivaslab data
external_files_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/imp/pgen/"
analysis_name = "rivaslab_ukbb_imputed_30k_rand_controls_sex_age"
covariates_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phe_qc/ukb24983_GWAS_covar.phe"
subject_ids_path = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/sqc/population_stratification/ukb24983_white_british.phe"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/"
mega_data_bim_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data.bim"

# # Read in sex, age, and PCs
# cov_data = read.delim(covariates_path)
# rownames(cov_data) = cov_data[,1]
# print("covariates loaded, dim: ")
# print(dim(cov_data))

chrs = 1:22
all_files = list.files(external_files_path)
if(analysis_name != ""){
  out_path = paste(out_path,analysis_name,"/",sep="")
  system(paste("mkdir",out_path))
}

script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

mega_bim = read.table(mega_data_bim_file,stringsAsFactors = F,header = F)
id_is_loc = grepl(":",mega_bim[,2])
mega_locations = paste(mega_bim[,1],mega_bim[,4],sep=":")
ids = mega_bim[,2]
names(mega_locations) = ids
corrected_locations = sapply(ids[id_is_loc],function(x)strsplit(x,split="-")[[1]][1])
mega_locations[id_is_loc] = corrected_locations

####################################################################################################
####################################################################################################
####################################################################################################
# Create frq files (useful for comparison with other mafs)
for (chr in chrs){
  curr_file = all_files[grepl(".bed$",all_files) & grepl(paste("chr",chr,"_",sep=""),all_files)]
  curr_file = gsub(pattern = ".bed$",replacement = "",curr_file)
  err_path = paste(out_path,"merge_geno.err",sep="")
  log_path = paste(out_path,"merge_geno.log",sep="")
  curr_cmd = paste("plink --bfile",paste(external_files_path,curr_file,sep=''),
                   "--keep",subject_ids_path,
                   "--freq --out",paste(out_path,"chr",chr,sep=""))
  curr_sh_file = "merge_geno.sh"
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}

# Comment this out: compare to our old mafs
# ukbb_hrc_data = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/FreqPlot-merged_control_geno-HRC.txt",
#                            stringsAsFactors = F)
# ukbb_frq = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq",
#                       stringsAsFactors = F,header=T)
# inds = which(abs(ukbb_hrc_data[,4]) > 0.05)
# ukbb_hrc_data[inds,]
# Comment this out: compare rivaslab freqs with our imputation
chr = 10
rivaslab_freq = read.table(paste(out_path,"chr",chr,".frq",sep=""),stringsAsFactors = F,header=T)
rownames(rivaslab_freq) = rivaslab_freq[,2]

# from direct imp folder
ashleylab_bim = read.table(paste("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/ukb_imp_chr",chr,"_v2.bim",sep=""),stringsAsFactors = F,
                           header=F)
ashleylab_freq = read.table(paste("/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/plink_small/freqs/chr",chr,".frq",sep=""),stringsAsFactors = F,
                           header=T)
# Our subset
ashleylab_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.bim",
                           stringsAsFactors = F,header=F)
ashleylab_freq = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq",
                            stringsAsFactors = F,header=T)

all(ashleylab_bim[,2]==ashleylab_freq[,2])
ashleylab_ids = paste(paste(ashleylab_bim[,1],ashleylab_bim[,4],sep=":"),ashleylab_bim[,5],ashleylab_bim[,6],sep="_")
rownames(ashleylab_freq) = ashleylab_ids

inds = intersect(ashleylab_ids,rivaslab_freq[,2])
r_mafs = rivaslab_freq[inds,"MAF"]
a_mafs = ashleylab_freq[inds,"MAF"]
cor(r_mafs,a_mafs)
diffs = r_mafs-a_mafs
table(diffs > 0.01)/sum(a_mafs | r_mafs > 0.01)
# # chr 2 output:
# > cor(r_mafs,a_mafs)
# [1] 0.9985048
# > diffs = r_mafs-a_mafs
# > table(diffs > 0.01)/length(diffs)
# FALSE       TRUE 
# 0.94054091 0.05945909 

# check one of our results: "22:35884154_G_A"
is.element(set=ashleylab_ids,"22:35884154_G_A")
ashleylab_freq["22:35884154_A_G",]
rivaslab_freq["22:35884154_G_A",]

# Check: rs143940620, 10:18054131_C_T 
# Or 10:105815241_A_C
ashleylab_freq["10:18054131_T_C",]
rivaslab_freq["10:18054131_T_C",]
# > ashleylab_freq["10:18054131_T_C",]
# CHR         SNP A1 A2    MAF NCHROBS
# 10:18054131_T_C  10 rs143940620  T  C 0.1033   26508
# > rivaslab_freq["10:18054131_T_C",]
# CHR             SNP A1 A2    MAF NCHROBS
# 10:18054131_T_C  10 10:18054131_T_C  T  C 0.1094  444406
# > 

# rsids: rs1747677;rs10748861
ashleylab_freq["10:105815241_A_C",]
rivaslab_freq["10:105815241_A_C",]
# > ashleylab_freq["10:105815241_A_C",]
# CHR       SNP A1 A2    MAF NCHROBS
# 10:105815241_A_C  10 rs1747677  A  C 0.3283   22314
# > rivaslab_freq["10:105815241_A_C",]
# CHR              SNP A1 A2    MAF NCHROBS
# 10:105815241_A_C  10 10:105815241_A_C  A  C 0.3081  375630
# cd /oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls
# less bfile1.bim | grep 105815241
# 10	rs1747677	122.922	105815241	A	C
# less bfile1.hwe | grep rs1747677
# 10  rs1747677  ALL(NP)    A    C   48/414/883   0.3078   0.3073    1
# less bfile1.frq | grep rs1747677
# 10    rs1747677    A    C       0.1896     2690
# Check freqs in recalling with genepool
# cd /oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall
# less plink_test_mega_consortium_data.frq | grep rs1747677
# rs1747677    A    C       0.3095     6928
# cd /oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/
# less FreqPlot-merged_mega_data_autosomal-HRC.txt | grep rs1747677
# rs1747677	0.340006	0.1901	0.149906	1

####################################################################################################
####################################################################################################
####################################################################################################
# Create reduced data files: ukbb whites and mega locations
for (chr in chrs){
  curr_file = all_files[grepl(".bed$",all_files) & grepl(paste("chr",chr,"_",sep=""),all_files)]
  curr_file = gsub(pattern = ".bed$",replacement = "",curr_file)
  
  # curr_bim  = read.table(paste(external_files_path,curr_file,".bim",sep=""),stringsAsFactors = F,header = F)
  # curr_locations = paste(curr_bim[,1],curr_bim[,4],sep=":")
  # to_rem = !is.element(curr_locations,set=mega_locations)
  # print(paste("looking at chromosome",chr,"number of locations shared with mega:",sum(!to_rem)))
  # print(paste("out of a total of mega snps:",sum(grepl(paste(chr,":",sep=""),mega_locations))))
  # curr_excluded = curr_bim[to_rem,2]
  curr_excluded_file = paste(out_path,"excluded_chr",chr,".txt",sep="")
  # write.table(t(t(curr_excluded)),quote = F,sep="",row.names = F,col.names = F,file=curr_excluded_file)
  
  err_path = paste(out_path,"merge_geno.err",sep="")
  log_path = paste(out_path,"merge_geno.log",sep="")
  curr_cmd = paste("plink --bfile",paste(external_files_path,curr_file,sep=''),
                   "--keep",subject_ids_path,
                   "--exclude",curr_excluded_file,
                   "--freq --make-bed --out",paste(out_path,"mega_snps_chr",chr,sep=""))
  curr_sh_file = "merge_geno.sh"
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}

