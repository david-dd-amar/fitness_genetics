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


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Locally: get the metadata and compare to Malene's metadata
setwd("/Users/David/Desktop/elite/metadata/")
load('idats_metadata.RData')
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
idat_paths = unname(sapply(unique(idat_files),function(x){
  arr=strsplit(x,split='/')[[1]];
  n = length(arr)
  paste(arr[1:(n-2)],collapse='/')
}))
unique(idat_paths)
table(idat_paths)/2

# Read metadata from different sources
mdata_elite = read.xlsx2('Elite_Cooper_metadata.xlsx',1)
mdata_cooper = read.xlsx2('Elite_Cooper_metadata.xlsx',2)
try({
	mdata_elite = read.xlsx2('fitness_genetics_ashleylab_metadata.xlsx',sheetName="ELITE_may18")
	mdata_cooper = read.xlsx2('fitness_genetics_ashleylab_metadata.xlsx',sheetName = "Cooper_may18")
})

# Look at Stanford3k
mdata_st3k = read.xlsx2("stanford3k_metadata2.xlsx",1)
try({
	mdata_st3k = read.delim("metadata/stanford3k_metadata2.txt",sep=" ")
})

# Same data from Kirstie
mdata_kirs = read.xlsx2("Stanford_Ashley_MEGAv2_n3484_DNAReport_Kirstie_dw.xlsx",1)
try({mdata_kirs = read.delim("metadata/Stanford_Ashley_MEGAv2_n3484_DNAReport_Kirstie_dw.tsv")})

# # We now extend the batch mdata files, we need it in this format for R packages
# samps_tmp = dnaid2sampleid[idats_dna_id]
# stanford3k_annot = unique(cbind(samps_tmp,idat_barcodes,idat_locs))
# colnames(stanford3k_annot) = colnames(mdata_batch1)
# rownames(stanford3k_annot) = NULL
# all_mdata = rbind(mdata_batch1,mdata_batch2)
# all_mdata = rbind(all_mdata,stanford3k_annot)
# rownames(all_mdata) = NULL
# write.csv(all_mdata,file="merged_metadata_file_stanford3k_elite_cooper.csv")

# # alternative: look at the dna ids directly and compare
# sample_ids_from_dna_ids = unlist(sapply(idats_dna_id,function(x){arr=strsplit(x,split='_')[[1]];arr[length(arr)]}))
# table(dnaid2sampleid[names(sample_ids_from_dna_ids)] == sample_ids_from_dna_ids)
# diff_dna_ids = which(dnaid2sampleid[names(sample_ids_from_dna_ids)] != sample_ids_from_dna_ids)
# sample_ids_from_dna_ids[names(diff_dna_ids)[200]]
# dnaid2sampleid[names(diff_dna_ids)[200]]

################# New script mapping sample ids to their idat files ##################
library(xlsx)
reverse_names<-function(x,remove_wo_names=T){
  v = names(x)
  names(v)=x
  if(remove_wo_names){
    to_rem = x==""|is.na(x)
    v = v[!to_rem]
  }
  return(v)
}
add_based_on_new_ids<-function(x,y,l1,l2){
  y = y[!is.element(y,set=rownames(x))]
  newx = cbind(l1[names(y)],l2[names(y)])
  rownames(newx) = unname(y)
  return(rbind(x,newx))
}
remap_rownames_using_alt_names<-function(x,y,alt){
  rows = rep(NA,length(x))
  names(rows) = x
  in_rownames = is.element(x,set=y)
  rows[in_rownames] = x[in_rownames]
  rows[!in_rownames] = alt[x[!in_rownames]]
  return(rows)
}

# Read the sample ids, including the alternate ids 
all_sample_data = read.xlsx2('fitness_genetics_ashleylab_metadata.xlsx',1,stringsAsFactors = F)
all_sample_data[all_sample_data$Sample_ID=="",1]=all_sample_data[all_sample_data$Sample_ID=="",]$Clinical_ID
sample_id_to_alternate_id = all_sample_data$alt_sample_id
names(sample_id_to_alternate_id) = all_sample_data$Sample_ID
sample_id_to_alternate_id = sample_id_to_alternate_id[sample_id_to_alternate_id!=""]
sample_id2dna_id_mainfile = all_sample_data$DNA_ID
names(sample_id2dna_id_mainfile) = all_sample_data$Sample_ID
clinical_id2sample_id = all_sample_data$Clinical_ID
names(clinical_id2sample_id) = all_sample_data$Sample_ID
clinical_id2sample_id = clinical_id2sample_id[clinical_id2sample_id!=""]
alternate_id2sample_id = c(alternate_id2sample_id,clinical_id2sample_id)
length(alternate_id2sample_id)
alternate_id2sample_id = reverse_names(sample_id_to_alternate_id)

# use the stanford3k file to map dna ids to sample ids
dnaid2sampleid_st3k = as.character(mdata_st3k$FID)
names(dnaid2sampleid_st3k) = mdata_st3k$IID

# This is the main table we create here. The addition rules below are ordered.
sample_id2idat_file = c()
unmapped_observed_idats = c()
unmapped_observed_dnaids = c()
# Rule 1: use the reports that use sample ids
mdata_batch1 = as.matrix(read.csv('1092_samples_report_1.csv',skip = 15,header = T)[,c(1,4,5)])
mdata_batch2 = as.matrix(read.csv('485_samples_report_2.csv',skip = 15,header = T)[,c(1,4,5)])
rownames(mdata_batch1) = as.character(mdata_batch1[,1])
rownames(mdata_batch2) = as.character(mdata_batch2[,1])
mdata_reports = rbind(mdata_batch1[,2:3],mdata_batch2[,2:3])
curr_rows_in_mainfile = remap_rownames_using_alt_names(rownames(mdata_reports),
                            all_sample_data$Sample_ID,alternate_id2sample_id)
curr_newrows = mdata_reports[!is.na(curr_rows_in_mainfile),]
rownames(curr_newrows) = curr_rows_in_mainfile[!is.na(curr_rows_in_mainfile)]
sample_id2idat_file = rbind(sample_id2idat_file,curr_newrows)
unmapped_observed_idats = rbind(unmapped_observed_idats,mdata_reports[is.na(curr_rows_in_mainfile),])

# Rule 2: use the DNA IDs in the mainfile to map to the idat files
dnaid2sampleid_mainfile = reverse_names(sample_id2dna_id_mainfile)
mainfile_shared_dnaids_with_idats = intersect(names(dnaid2sampleid_mainfile),names(idat_barcodes))
curr_dnaid2sample_id = dnaid2sampleid_mainfile[mainfile_shared_dnaids_with_idats]
# as of may 2018 this is empty:
unmapped_but_in_mainfile = names(which(""==curr_dnaid2sample_id | is.na(curr_dnaid2sample_id)))
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,curr_dnaid2sample_id,idat_barcodes,idat_locs)
dim(sample_id2idat_file)
table(is.element(rownames(sample_id2idat_file),set=all_sample_data$Sample_ID))
table(is.element(rownames(sample_id2idat_file),set=names(alternate_id2sample_id)))

# Rule 3: use stanford 3k DNA IDs in the mainfile to map to the idat files
st3k_shared_dnaids_with_idats = intersect(names(dnaid2sampleid_st3k),names(idat_barcodes))
curr_dnaid2sample_id = dnaid2sampleid_st3k[st3k_shared_dnaids_with_idats]
table(is.element(curr_dnaid2sample_id,set=names(alternate_id2sample_id)))
table(is.element(curr_dnaid2sample_id,set=all_sample_data$Sample_ID))
# As of may 2018: there are no alternative ids here
curr_dnaid2sample_id = curr_dnaid2sample_id[is.element(curr_dnaid2sample_id,set=all_sample_data$Sample_ID)]
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,curr_dnaid2sample_id,idat_barcodes,idat_locs)
dim(sample_id2idat_file)
table(is.element(rownames(sample_id2idat_file),set=all_sample_data$Sample_ID))
table(is.element(rownames(sample_id2idat_file),set=names(sample_id_to_alternate_id)))

# Rule 4: use the dna ids in the idat file names
sample_ids_from_dna_ids = unlist(sapply(idats_dna_id,function(x){arr=strsplit(x,split='_')[[1]];arr[length(arr)]}))
sample_ids_from_dna_ids = sample_ids_from_dna_ids[is.element(sample_ids_from_dna_ids,
    set=union(all_sample_data$Sample_ID,names(alternate_id2sample_id)))]
curr_alt_ids = is.element(sample_ids_from_dna_ids,set=names(alternate_id2sample_id))
table(curr_alt_ids)
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,sample_ids_from_dna_ids,idat_barcodes,idat_locs)
dim(sample_id2idat_file)

# Get stats by cohorts
curr_rows = is.element(all_sample_data$Sample_ID,set=rownames(sample_id2idat_file))
table(all_sample_data$Cohort[curr_rows])
table(table(sample_id2idat_file[,1]))
length(setdiff(unique(idat_barcodes),unique(sample_id2idat_file[,1])))
length(intersect(unique(idat_barcodes),unique(sample_id2idat_file[,1])))



