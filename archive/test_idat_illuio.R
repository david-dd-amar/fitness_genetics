# This script goes over idat files and read them using illuminaio
RLIBS_path = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/rpackages/"
WD_path = "/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/"
# within the WD
data_paths = c(
  "raw/stanford3k/idats/",
  "raw/idats_may_2018",
  "raw/elite_batch_2019",
  "raw/elite_batch_2020"
)

setwd(WD_path)

# Older R versions
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("illuminaio",lib="rpackages")
# # biocLite("crlmm",lib=RLIBS_path)
# library("oligoClasses",lib=RLIBS_path)
# library("preprocessCore",lib=RLIBS_path)
# library("illuminaio",lib=RLIBS_path)
# library("crlmm",lib=RLIBS_path)

# On R 3.6.0
# Preinstalled in R, install if required
cran_packages = c(
  "bitops","RCurl","rJava","openssl","bit64","RSQLite"
)
for(pkg in cran_packages){
  library(pkg,character.only = T)
}

# When installing packages below, tell R to update the outdated packages
# (R will ask you whether to update or not)
# install.packages("BiocManager",lib=RLIBS_path,dependencies = T)
bioc_packages = c("GenomicRanges","GenomeInfoDb",
                  "BiocGenerics","zlibbioc",
                  "affyio","Biobase","IRanges","oligoClasses",
                  "preprocessCore","illuminaio")
# for(pkg in bioc_packages){
#   BiocManager::install(pkg,lib = RLIBS_path,dependencies=T)
# }

for(pkg in bioc_packages){
  library(pkg,lib.loc = RLIBS_path,character.only = T)
}

data_path_idats = list()
for(p in data_paths){
  cmd = paste("find",p,"|grep \"idat$\" > tmp.txt")
  print(cmd)
  system(cmd)
  data_path_idats[[p]] = readLines("tmp.txt")
  system("rm tmp.txt")
}
sapply(data_path_idats,length)/2

illuminaio_extract_idat_metadata<-function(path){
  obj = illuminaio::readIDAT(path)
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
	print(path)
	return(list(
	  path=path, unknowns = ukns,
		versionNum = version, date=date_scanned,
		nSNPs=nSNPs,is_red=is_red,
		chip_type=chipType
	))
}

try_idat_read_several_times<-function(path,tries = 10){
  for(j in 1:tries){
    res = illuminaio_extract_idat_metadata(path)
    if(length(res)==7){return(res)}
  }
  return(NULL)
}
# illuminaio_extract_idat_metadata(data_path_idats[[1]][[1]])

# Read all idats and keep the metadata
idats_metadata = list()
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
for(nn1 in names(data_path_idats)){
  print(nn1)
  idats_metadata[[nn1]] = mclapply(data_path_idats[[nn1]],illuminaio_extract_idat_metadata,mc.cores = 32)
  names(idats_metadata[[nn1]]) = data_path_idats[[nn1]]
  save(idats_metadata,file="idats_metadata.RData")
}

sapply(idats_metadata,length) == sapply(data_path_idats,length)
lengths = sapply(idats_metadata,function(x)sapply(x,length))
sapply(lengths,table)
all(sapply(idats_metadata[[1]],function(x)x$path) == names(idats_metadata[[1]]))

# Relevant only if there are NULL objects that were read
# fail_examples = which(lengths[[4]]==0)[1:5]
# data_path_idats[[4]][fail_examples]

# Compare the may 2018 new batch to the old one from 2017
library(tools)
path2017 = '/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/raw/EA171116.01.0-02.1-delivered'
path2018 = '/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/raw/idats_may_2018/iData_Euan_Ashley'
files2017 = sort(list.dirs(path2017,full.names=T,recursive=F))
names(files2017) = sort(list.dirs(path2017,full.names=F,recursive=F))
files2018 = sort(list.dirs(path2018,full.names=T,recursive=F))
names(files2018) = sort(list.dirs(path2018,full.names=F,recursive=F))
print(length(intersect(names(files2017),names(files2018))))
print(length(setdiff(names(files2017),names(files2018))))
print(length(setdiff(names(files2018),names(files2017))))
print(length(union(names(files2017),names(files2018))))
for(nn in intersect(names(files2017),names(files2018))){
	files1 = list.files(files2017[nn],full.names=T)
	names(files1) = list.files(files2017[nn],full.names=F)
	files2 = list.files(files2018[nn],full.names=T)
	names(files2) = list.files(files2018[nn],full.names=F)
	print(all(names(files1)==names(files2)))
	comp = md5sum(files1)==md5sum(files2)
	print(table(comp))
}

# Create a file with paths, sample ids, and batches
load("idats_metadata.RData")
df = c()
for(nn in names(idats_metadata)){
  print(nn)
  x = sapply(idats_metadata[[nn]],unlist)
  x = t(x)
  rownames(x) = NULL
  print(x[1:5,2:8])
  df = rbind(df,x)
}


