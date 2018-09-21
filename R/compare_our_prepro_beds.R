# In this script we compare two bed files that shared subjects.
# The main usage is to compare our preprocessing from genomestudio to
# the reports we have from Illumina.
# Assumption: beds where created using the same strand annotation.
# Nevertheless, we go over the bim and take identical snps only.

script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Locations of old plink files from Illumina's analysis:
# 1. The consortium chip data: there are many options, examples below
#     /scratch/groups/euan/projects/stanford3k/plink/qc/mega_n3484_filter
#     /scratch/groups/euan/projects/stanford3k/plink/binary/Stanford_Ashley_MEGAv2_n3484_gtReport_File-1_top
# 2. The mega global array (other batches):
#     /oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/batch1_1092/sample_report_1092
# We have either plus strand or TOP/BOT strand called plink data

# Keep in mind! We use plink only to find identical or nearly identical samples.
# We therefore do not prune the snps and we also sample the data for running time. 
# The actual ibd scores are not usefule per se. They are only useful for detecting
# "twins" or repeated samples by getting an almost perfect fit, i.e., by "distance" > 0.999
# or PI_HAT ~ 1.

bfile1 = "/scratch/groups/euan/projects/stanford3k/plink/qc/mega_n3484_filter"
# bfile1 = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/batch1_1092/sample_report_1092/samples_report_1092.95.snps.clean"

# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/raw" # our data
# bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering" # our data
bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/plink_test_mega_consortium_data/test_mega_consortium"

out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/qc/bedchecks/mega_comparisons/"
metadata_file = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"

bim_data1 = read.table(paste(bfile1,".bim",sep=""),stringsAsFactors = F,header = F)
rownames(bim_data1) = as.character(bim_data1[,2])

bim_data2 = read.table(paste(bfile2,".bim",sep=""),stringsAsFactors = F,header = F)
rownames(bim_data2) = as.character(bim_data2[,2])

shared_snps = intersect(rownames(bim_data1),rownames(bim_data2))
print(paste("number of snps in bfile 1:",nrow(bim_data1)))
print(paste("number of snps in bfile 2:",nrow(bim_data2)))
print(paste("number of shared snp ids:",length(shared_snps)))

# Check if snps are actually the same
check_snps = bim_data1[shared_snps,5]==bim_data2[shared_snps,5] &
  bim_data1[shared_snps,6]==bim_data2[shared_snps,6]
print(paste("number of shared snp ids with the same snp info:",sum(check_snps)))
names(check_snps) = shared_snps
shared_snps = names(which(check_snps))

# Exclude JHU snps
shared_snps = shared_snps[!grepl("^JHU",shared_snps)]

# Go over the mapping of the different ids
d2 = read.delim(metadata_file,stringsAsFactors = F)
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID;altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids; names(altsamp_id) = d2_ids
clin_id = d2$Clinical_ID
names(clin_id) = d2_ids


# Check subject overlap
fam_data1 = read.table(paste(bfile1,".fam",sep=""),stringsAsFactors = F,header = F)
sample_ids_bfile1 = sapply(fam_data1[,2],function(x){arr = strsplit(x,split="_")[[1]];arr[length(arr)]})
rownames(fam_data1) = as.character(fam_data1[,2])
fam_data2 = read.table(paste(bfile2,".fam",sep=""),stringsAsFactors = F,header = F)
rownames(fam_data2) = as.character(fam_data2[,2])

# fam_data2_sampid = samp_id[rownames(fam_data2)]
# shared_subjects = intersect(sample_ids_bfile1,fam_data2_sampid)
# print(paste("number of subjects in bfile 1:",nrow(fam_data1)))
# print(paste("number of subjects in bfile 2:",nrow(fam_data2)))
# print(paste("number of shared subject ids:",length(shared_subjects)))

# A better way to get subject overlap: get the ids directly from the idat files
load('/oak/stanford/groups/euan/projects/fitness_genetics/metadata/idats_metadata.RData')
idats_metadata_table = t(sapply(idats_metadata,function(x)x))
idats_dna_id = unname(idats_metadata_table[,5])
idats_dna_id = gsub(idats_dna_id,pattern=" ",replace="")
idats_dna_id = gsub(idats_dna_id,pattern="-DNA",replace="_")
idats_dna_id = gsub(idats_dna_id,pattern="__",replace="_")
idat_files = rownames(idats_metadata_table)
idat_barcodes = unname(sapply(idat_files,function(x){arr=strsplit(x,split='/')[[1]];arr[length(arr)-1]}))
idat_locs = unname(idats_metadata_table[,3])
names(idat_barcodes) = idats_dna_id
names(idat_locs) = idats_dna_id
our_ids = paste(idat_barcodes,idat_locs,sep="_")
names(our_ids) = idats_dna_id
names(idats_dna_id) = our_ids

fam2_ids_to_fam1 = c()
for(id in fam_data2[,2]){
	dnaid = idats_dna_id[id]
	arr = strsplit(dnaid,split="_")[[1]]
	regex1 = arr[1]
	regex2 = paste(arr[-1],collapse="_")
	inds = grepl(regex1,fam_data1[,2]) & grepl(regex2,fam_data1[,2])
	if(sum(inds)!=1){next}
	if(sum(inds)>0){
		fam2_ids_to_fam1[id] = fam_data1[inds,2][1]
	}
}
print(paste("number of subjects in bfile 1:",nrow(fam_data1)))
print(paste("number of subjects in bfile 2:",nrow(fam_data2)))
print(paste("number of shared subject ids:",length(fam2_ids_to_fam1)))

N1 = 22
N2 = 100000 
R = 1
for(j in 1:R){
  curr_subjects = sample(fam2_ids_to_fam1)[1:N1]
  curr_snps = sample(shared_snps)[1:N2]
  shared_subjects_bfile1 = curr_subjects
  print(table(shared_subjects_bfile1))
  shared_subjects_bfile2 = names(curr_subjects)
  print(table(shared_subjects_bfile2))
  currp = paste(out_path,"tmp",j,"/",sep="")
  system(paste("mkdir",currp))
  write.table(fam_data1[shared_subjects_bfile1,1:2],file=paste(currp,"subj1.txt",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  write.table(fam_data2[shared_subjects_bfile2,1:2],file=paste(currp,"subj2.txt",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  write.table(t(t(curr_snps)),file=paste(currp,"snps.txt",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  
  err_path = paste(currp,"tmp_bfile1_creation.err",sep="")
  log_path = paste(currp,"tmp_bfile1_creation.log",sep="")
  curr_cmd = paste("plink --bfile",bfile1,
                   "--keep",paste(currp,"subj1.txt",sep=''),
                   "--extract",paste(currp,"snps.txt",sep=""),
                   "--recode --make-bed --out",paste(currp,"bfile1",sep=''))
  curr_sh_file = "tmp_bfile1_creation.sh"
  print_sh_file(paste(currp,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(currp,curr_sh_file,sep='')))
  
  err_path = paste(currp,"tmp_bfile2_creation.err",sep="")
  log_path = paste(currp,"tmp_bfile2_creation.log",sep="")
  curr_cmd = paste("plink --bfile",bfile2,
                   "--keep",paste(currp,"subj2.txt",sep=''),
                   "--extract",paste(currp,"snps.txt",sep=""),
                   "--recode --make-bed --out",paste(currp,"bfile2",sep=''))
  curr_sh_file = "tmp_bfile2_creation.sh"
  print_sh_file(paste(currp,curr_sh_file,sep=''),
                get_sh_default_prefix(err_path,log_path),curr_cmd)
  system(paste("sbatch",paste(currp,curr_sh_file,sep='')))
  
  wait_for_job()
  
  # merge and run ibd
  err_path = paste(currp,"merge_plink.err",sep="")
  log_path = paste(currp,"merge_plink.log",sep="")
  curr_cmd = paste("plink --bfile",paste(currp,"bfile1",sep=''),
                   "--bmerge",paste(currp,"bfile2",sep=''),
                   "--distance square ibs --genome",
                   "--make-bed --out",paste(currp,"merged_data",sep=''))
  curr_sh_file = "merge_plink.sh"
  print_sh_file(paste(currp,curr_sh_file,sep=''),
                get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=1,mem_size=4000),curr_cmd)
  system(paste("sbatch",paste(currp,curr_sh_file,sep='')))
}
wait_for_job()

# Go over the results and compare the matching samples
for(j in 1:R){
  currp = paste(out_path,"tmp",j,"/",sep="")
  fam_info = read.table(paste(currp,"merged_data.fam",sep = ""),stringsAsFactors = F)
  N = nrow(fam_info)
  ids1 = 1:(N/2)
  ids2 = ((N/2)+1):N
  currdist = as.matrix(read.table(paste(currp,"merged_data.mibs",sep="")))
  
  # read the peds directly
  currp = paste(out_path,"tmp",j,"/",sep="")
  ped1 = read.table(paste(currp,"bfile1.ped",sep=""),stringsAsFactors=F)
  
  paired_dists = c()
  for(i in ids1){
	id1 = fam_info[i,2]
	id2 = fam2_ids_to_fam1[id1]
    i2 = which(fam_info[,2] == id2)
    paired_dists = rbind(paired_dists,c(i,i2,currdist[i,i2]))
  }
  print(mean(paired_dists[,3]))
}












