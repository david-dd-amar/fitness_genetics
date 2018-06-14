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
for(p in data_paths[4]){
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

# Compare the may 2018 new batch to the old one from 2017
library(tools)
path2017 = '/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/raw/EA171116.01.0-02.1-delivered'
path2018 = '/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/raw/idats_may_2018/iData_Euan_Ashley'
files2017 = sort(list.dirs(path2017,full.names=T,recursive=F))
names(files2017) = sort(list.dirs(path2017,full.names=F,recursive=F))
files2018 = sort(list.dirs(path2018,full.names=T,recursive=F))
names(files2018) = sort(list.dirs(path2018,full.names=F,recursive=F))
for(nn in intersect(names(files2017),names(files2018))){
	files1 = list.files(files2017[nn],full.names=T)
	names(files1) = list.files(files2017[nn],full.names=F)
	files2 = list.files(files2018[nn],full.names=T)
	names(files2) = list.files(files2018[nn],full.names=F)
	print(all(names(files1)==names(files2)))
	comp = md5sum(files1)==md5sum(files2)
	print(table(comp))
}


