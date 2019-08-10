# In this script we take two paths with bed files and merge each chromosome
# Assumption, check_bim.pl was run on each bed.
# As for the input, note that the analysis is not symmetric. Whenever possible, we convert the ids
# in path2 to match those in path1 For example if we use ukbb for path1 and Illumina ids for 
# path2 then non-standard ids such as (examXXX) will be mapped to ids from ukbb.
# Merge is done using plink
# Output: a merged bed file for each chromosome

# March 2019: MEGA with GP, Impute2, 1000G vs. our imputed UKBB
path1 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/ukbb_20k_imp/impute2_1000gRef_out/check_bim_res/"
path2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/impute2_1000gRef_out/check_bim_res/"
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/with_our_imp_ukbb/"

# define paths and libs
try(system(paste("mkdir",out_path)))
check_bim_info = T
qctool_path = "/home/users/davidama/apps/qctool_v2/build/release/qctool_v2.0.1"
script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)
library(data.table,lib.loc = "~/R/packages")

maf_for_pca = 0.01
exclude_pali = T
chrs = paste("chr",c(1:22),sep="")
is_snp_paly<-function(x){
  return(x=="AT" || x=="TA" || x=="CG" || x=="GC")
}

for(chr in chrs){
  print(paste("analyzing chr:",chr))
  # define the beds
  bfile1 = paste(path1,chr,sep="")
  bfile2 = paste(path2,chr,sep="")
  
  # Compare the bim files
  # (1) Check SNP intersect, locations, and which snps must be flipped before we analyze
  bim_data1 = fread(paste(bfile1,".bim",sep=""),stringsAsFactors = F,header = F,data.table = F)
  bim_data2 = fread(paste(bfile2,".bim",sep=""),stringsAsFactors = F,header = F,data.table = F)
  rownames(bim_data2) = bim_data2[,2]
  rownames(bim_data1) = bim_data1[,2]
  
  # Get the shared SNPs
  shared_snps = intersect(bim_data2[,2],bim_data1[,2])
  curr_loc_intersect = intersect(bim_data1[,4],bim_data2[,4])
  alt_ids1 = apply(bim_data1[,4:6],1,paste,collapse=":")
  alt_ids2 = apply(bim_data2[,4:6],1,paste,collapse=":")
  alt_ids2_sanity = apply(bim_data2[,c(4,6,5)],1,paste,collapse=";")
  print(paste("intersection:",length(intersect(alt_ids1,alt_ids2))))
  print(paste("number of shared locations:",length(curr_loc_intersect)))
  # the following should be zero because of our assumption that
  # both datasets went through check_bim.pl
  print(paste("intersection sanity check (should be zero):",length(intersect(alt_ids1,alt_ids2_sanity))))
  alt_intersect = intersect(alt_ids1,alt_ids2)
  rownames(bim_data1) = alt_ids1
  rownames(bim_data2) = alt_ids2
  # another sanity check: all id-based shared snps should be shared in the 
  # new ids as well
  print(paste("sanity check2:",all(is.element(shared_snps,bim_data1[alt_intersect,2]))))
  bim_data1 = bim_data1[alt_intersect,];bim_data2 = bim_data2[alt_intersect,]
  curr_shared_rsids = cbind(bim_data1[,2],bim_data2[,2])
  rownames(curr_shared_rsids) = alt_intersect
  
  # Order the SNPs by location: prevents issues with PLINK
  ord = order(as.numeric(bim_data2[,4]))
  curr_shared_rsids = curr_shared_rsids[ord,]
  
  # In case we want to remove the palindromic SNPs
  if(exclude_pali){
    alleles = paste(bim_data2[,5],bim_data2[,6],sep="")
    pali_snps = sapply(alleles,is_snp_paly)
    print(paste("percent of pali snps:",sum(pali_snps)/length(pali_snps)))
    curr_shared_rsids = curr_shared_rsids[!pali_snps,]
  }
  
  # Print force-allele files - we are going to load the data to plink and change
  # some of the bim info. Here we make sure that we keep the input alleles before we 
  # extract and work with the shared snps.
  force_allele1 = bim_data1[rownames(curr_shared_rsids),c(2,5)]
  write.table(force_allele1,file=paste(out_path,chr,"_force_allele1.txt",sep=""),
              row.names = F,quote=F,col.names = F,sep="\t")
  force_allele1_cmd = paste("--reference-allele",paste(out_path,chr,"_force_allele1.txt",sep=""))
  force_allele2 = bim_data2[rownames(curr_shared_rsids),c(2,5)]
  write.table(force_allele2,file=paste(out_path,chr,"_force_allele2.txt",sep=""),
              row.names = F,quote=F,col.names = F,sep="\t")
  force_allele2_cmd = paste("--reference-allele",paste(out_path,chr,"_force_allele2.txt",sep=""))
  # sanity_check
  print(paste("sanity check3:",all(as.character(force_allele1[,2])==as.character(force_allele2[,2]))))
  
  extract_snps_using_plink(bfile2,curr_shared_rsids[,2],out_path,
                           paste(chr,"_file2_shared_snps",sep=""),
                           paste(chr,"_new_bed_2",sep=""),
                           ref_allale_line = force_allele2_cmd,
                           get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
  extract_snps_using_plink(bfile1,curr_shared_rsids[,1],out_path,
                           paste(chr,"_file1_shared_snps",sep=""),
                           paste(chr,"_new_bed_1",sep=""),
                           ref_allale_line = force_allele1_cmd,
                           get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000)
}
wait_for_job(waittime = 60)

total_num_snps = 0
for(chr in chrs){
  # Compare the new two beds
  new_bed_2_bim = fread(paste(out_path,chr,"_new_bed_2.bim",sep=""),
                        stringsAsFactors = F,data.table = F)
  rownames(new_bed_2_bim) = new_bed_2_bim[,2]
  new_bed_1_bim = fread(paste(out_path,chr,"_new_bed_1.bim",sep=""),
                        stringsAsFactors = F,data.table=F)
  rownames(new_bed_1_bim) = new_bed_1_bim[,2]
  print(paste("sanity checks, all should be true:"))
  print(all(new_bed_1_bim[,4]==new_bed_2_bim[,4]))
  print(!is.unsorted(new_bed_1_bim[,4]))
  print(all(as.character(new_bed_1_bim[,5])==as.character(new_bed_2_bim[,5])))
  print(all(as.character(new_bed_1_bim[,6])==as.character(new_bed_2_bim[,6])))
  diff_rows = which(as.character(new_bed_1_bim[,6])!=as.character(new_bed_2_bim[,6]))
  print(length(diff_rows)==0)
  
  # Print the new bim file for file 2
  new_names = new_bed_1_bim[,2]
  new_bed_2_bim[,2] = new_names
  rownames(new_bed_2_bim) = NULL
  write.table(new_bed_2_bim,
              file=paste(out_path,chr,"_new_bed_2_alt.bim",sep=""),sep="\t",
              row.names=F,col.names=F,quote=F)
  total_num_snps = total_num_snps+nrow(new_bed_2_bim)
  
  # Merge using PLINK
  err_path = paste(out_path,chr,".err",sep="")
  log_path = paste(out_path,chr,".log",sep="")
  curr_cmd = paste("plink --bed",paste(out_path,chr,"_new_bed_2.bed",sep=''),
                   "--bim",paste(out_path,chr,"_new_bed_2_alt.bim",sep=''),
                   "--fam", paste(out_path,chr,"_new_bed_2.fam",sep=''),
                   "--bmerge",paste(out_path,chr,"_new_bed_1",sep=''),
                   "--reference-allele",paste(out_path,chr,"_force_allele1.txt",sep=""),
                   "--threads 4",
                   "--make-bed --out",paste(out_path,chr,sep=''))
  curr_sh_file = paste(out_path,chr,".sh",sep="")
  print_sh_file(curr_sh_file,
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),
                curr_cmd)
  system(paste("sbatch",curr_sh_file))  
}


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# Run flipscan: assumption - we now have a set of merged bfiles, one per chromosome

fam1 = as.matrix(read.table(paste(bfile1,".fam",sep=""),stringsAsFactors = F))
fam2 = as.matrix(read.table(paste(bfile2,".fam",sep=""),stringsAsFactors = F))
fam1 = cbind(fam1,rep("1",nrow(fam1)))
fam2 = cbind(fam2,rep("2",nrow(fam2)))
rownames(fam1)=NULL;colnames(fam1)=NULL
rownames(fam2)=NULL;colnames(fam2)=NULL
phe = rbind(fam1,fam2)
phe = phe[,c(1:2,7)]
colnames(phe) = c("FID","IID","merged")
write.table(phe,file=paste(out_path,"flipscan_type.txt",sep=""),
            sep=" ",quote=F,row.names = F,col.names = T)

# Run flipscan
for(chr in chrs){
  err_path = paste(out_path,chr,"_flipscan.err",sep="")
  log_path = paste(out_path,chr,"_flipscan.log",sep="")
  curr_cmd = paste("plink --bfile",paste(out_path,chr,sep=''),
                   "--flip-scan --allow-no-sex",
                   "--pheno",paste(out_path,"flipscan_type.txt",sep=""),
                   "--pheno-name merged",
                   "--out",paste(out_path,chr,"_merged_data_flipscan",sep=''))
  curr_sh_file = paste(chr,"_mega_flipscan.sh",sep="")
  print_sh_file(paste(out_path,curr_sh_file,sep=''),
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=4,mem_size=32000),
                curr_cmd)
  system(paste("sbatch",paste(out_path,curr_sh_file,sep='')))
}
wait_for_job(120)

# Read the flipscan results, exclude snps
excluded_snps_flipscan = c()
for(chr in chrs){
  flipscan_res = readLines(paste(out_path,chr,"_merged_data_flipscan.flipscan",sep=''))
  arrs = strsplit(flipscan_res[-1],split="\\s+")
  names(arrs) = sapply(arrs,function(x)x[3])
  flipscan_failures = sapply(arrs,length) > 11
  flipscan_failures = sapply(arrs[flipscan_failures],function(x)x[3])
  bim = fread(paste(out_path,chr,".bim",sep=""),stringsAsFactors = F,
              header=F,data.table = F)
  snps_to_keep = setdiff(bim[,2],flipscan_failures)
  excluded_snps_flipscan = c(excluded_snps_flipscan,flipscan_failures)
  print(paste("Flipscan check, number of variants to remove:",length(flipscan_failures)))
  extract_snps_using_plink(paste(out_path,chr,sep=''), # bfile
                           snps_to_keep,out_path, # list of snps to keep, output path
                           paste(chr,"_final_snps_to_keep_after_flipscan",sep=""), # the name for the tmp txt file
                           paste(chr,sep=''), # output bfile
                           batch_script_func=get_sh_prefix_one_node_specify_cpu_and_mem,
                           Ncpu=4,mem_size=32000)
}
save(excluded_snps_flipscan,file=paste(out_path,"excluded_snps_flipscan.RData",sep=""))
wait_for_job(60)



