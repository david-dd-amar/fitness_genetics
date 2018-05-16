# This script takes all the paths form the project's README,
# searches for all idats and read them using illuminaio
README_path = "/oak/stanford/groups/euan/projects/fitness_genetics/README"
RLIBS_path = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/rpackages/"
WD_path = "/oak/stanford/groups/euan/projects/fitness_genetics/"

setwd(WD_path)
# source("https://bioconductor.org/biocLite.R")
# biocLite("illuminaio",lib="rpackages")
# biocLite("crlmm",lib=RLIBS_path)
library("oligoClasses",lib=RLIBS_path)
library("preprocessCore",lib=RLIBS_path)
library("illuminaio",lib=RLIBS_path)
library("crlmm",lib=RLIBS_path)

# Get paths from README
data_paths = c()
for(l in readLines(README_path)){
  if(!grepl("PATH=",l)){next}
  data_paths=c(data_paths,strsplit(l,split="PATH=")[[1]][2])
}

data_path_idats = list()
for(p in data_paths[-4]){
  cmd = paste("find",p,"|grep \"idat$\" > tmp.txt")
  print(cmd)
  system(cmd)
  data_path_idats[[p]] = readLines("tmp.txt")
  print(data_path_idats[[p]])
  system("rm tmp.txt")
}
data_path_idats = data_path_idats[sapply(data_path_idats,length)>0]

#idat_files = as.character(read.table("analysis/all_idat_files.txt")[,1])
#batches = sapply(idat_files,function(x)strsplit(x,split="_")[[1]][1])

illuminaio_extract_idat_metadata<-function(obj){
	rem = which(sapply(obj,length)>500)
	obj = obj[-rem]
	ukns = unname(unlist(obj$Unknowns))
	version = obj$versionNumber
	date_scanned = obj$RunInfo[1,1]
	is_red = obj$RedGreen
	barcode = obj$Barcode
	nSNPs = obj$nSNPs
	chipType = obj$ChipType
	#info_table = obj$RunInfo
	#field_table = obj$fields
	# learn our sample id 
	ind = which(grepl('-',ukns) & grepl("_",ukns))
	sample_id = ""
	if(length(ind)>0){
		sample_id = strsplit(ukns[ind[1]],split='_')[[1]][2]
	}
	return(c(
		sample_id = sample_id,
		ukns,versionNum=version,date=date_scanned,nSNPs=nSNPs,is_red=is_red,
		chip_type=chipType
	))
}

# Read all idats and keep the metadata
idats_metadata = list()
try({load("idats_metadata.RData")})
for(nn1 in names(data_path_idats)){
	for (nn2 in data_path_idats[[nn1]]){
		if(is.element(nn2,set=names(idats_metadata))){next}
		obj = readIDAT(nn2)
		mdata = illuminaio_extract_idat_metadata(obj)
		idats_metadata[[nn2]] = c(mdata,ashleylab_batch_name = nn1)
		save(idats_metadata,file="idats_metadata.RData")
		print(length(idats_metadata))
	}
}

# Locally: get the metadata and compare to Malene's metadata
setwd("/Users/David/Desktop/elite/")
library(xlsx)
load('idats_metadata.RData')
idats_metadata_table = t(sapply(idats_metadata,function(x)x))
idats_dna_id = unname(idats_metadata_table[,5])
idats_dna_id = gsub(idats_dna_id,pattern=" ",replace="")
idats_dna_id = gsub(idats_dna_id,pattern="-DNA",replace="_")
idats_dna_id = gsub(idats_dna_id,pattern="__",replace="_")
idat_files = rownames(idats_metadata_table)
idat_barcodes = unname(sapply(idat_files,function(x){arr=strsplit(x,split='/')[[1]];arr[length(arr)-1]}))
idat_locs = unname(idats_metadata_table[,3])
mdata_elite = read.xlsx2('Elite_Cooper_metadata.xlsx',1)
mdata_cooper = read.xlsx2('Elite_Cooper_metadata.xlsx',2)
idat_paths = unname(sapply(unique(idat_files),function(x){
  arr=strsplit(x,split='/')[[1]];
  n = length(arr)
  paste(arr[1:(n-2)],collapse='/')
  }))
unique(idat_paths)
table(idat_paths)/2

# check the samples by id names
mdata_dna_ids = union(mdata_cooper$DNA_ID,mdata_elite$DNA_ID)
length(intersect(mdata_dna_ids,idats_dna_id))

# Look at Stanford3k
mdata_st3k = read.xlsx2("stanford3k_metadata2.xlsx",1)
length(intersect(mdata_st3k$IID,idats_dna_id))
length(unique(idats_dna_id))
length(unique(mdata_st3k$IID))
# use the stanford3k file to map dna ids to sample ids
dnaid2sampleid = as.character(mdata_st3k$FID)
names(dnaid2sampleid) = mdata_st3k$IID

# Same data from Kirstie
mdata_kirs = read.xlsx2("Stanford_Ashley_MEGAv2_n3484_DNAReport_Kirstie_dw.xlsx",1)
# idats internal ids vs the metadata
length(intersect(mdata_kirs$DNA_ID,idats_dna_id))
# Cooper?
cooper_sample_ids = mdata_cooper$Clinical_ID
length(intersect(mdata_kirs$Notes,cooper_sample_ids))

# Add data from the new batches
mdata_batch1 = read.csv('metadata/1092_samples_report_1.csv',skip = 15,header = T)[,c(1,4,5)]
mdata_batch2 = read.csv('metadata/485_samples_report_2.csv',skip = 15,header = T)[,c(1,4,5)]
new_batches_barcodes = union(mdata_batch1[,2],mdata_batch2[,2])
table(is.element(idat_barcodes,set=new_batches_barcodes))/2
table(is.element(idat_barcodes,set=mdata_batch1[,2]))/2

# We now extend the batch mdata files, we need it in this format for R packages
samps_tmp = dnaid2sampleid[idats_dna_id]
stanford3k_annot = unique(cbind(samps_tmp,idat_barcodes,idat_locs))
colnames(stanford3k_annot) = colnames(mdata_batch1)
all_mdata = rbind(mdata_batch1,mdata_batch2)
all_mdata = rbind(all_mdata,stanford3k_annot)
write.csv(all_mdata,file="merged_metadata_file_stanford3k_elite_cooper.csv")

# Test data locally
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("crlmm")
library(crlmm)
test_barcode = "200206160028"
curr_mdata = as.matrix(all_mdata[all_mdata[,2]==test_barcode,])
curr_mdata = data.frame(curr_mdata,stringsAsFactors = F)
write.csv(curr_mdata,file='idats/manifest.csv')
idata = readIdatFiles(curr_mdata,path="idats/")
class(idata)
idata[1,1]
gdata = genotype.Illumina(curr_mdata,path='idats')

library(devtools)
## allow R to look for pacakges in both CRAN and Bioconductor
setRepositories(ind = 1:2)
## install from Github source
install_github("andrewparkermorgan/argyle")
