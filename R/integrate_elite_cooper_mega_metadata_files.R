
# In this script we go over our lab's different annotation files.
# Our goal here is to map the DNA ids imprinted in the idat files to our
# sample ids. Note that we do not analyze the genepool samples in detail. 
# For this we have the script our_dataset_merge_with_gp_annotation, which
# takes the output of this script and merge it with the different genepool
# annotation sources. 
#
# Additional comments are given in line.

# Locally: get the metadata and compare to Malene's metadata
setwd("/Users/David/Desktop/elite/metadata/")
options(java.parameters = "-Xmx8000m")
library(xlsx)

# Load idat metadata
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
# # tests and stats
# idat_paths = unname(sapply(unique(idat_files),function(x){
#   arr=strsplit(x,split='/')[[1]];
#   n = length(arr)
#   paste(arr[1:(n-2)],collapse='/')
# }))
# unique(idat_paths)
# table(idat_paths)/2
# idats_dna_id[idat_barcodes=="200200980068" & idat_locs=="R02C01"]
# idats_dna_id[idat_barcodes=="201557580009" & idat_locs=="R02C01"]
# idats_dna_id[grepl("1391",idats_dna_id)]
# idats_dna_id[grepl("1390",idats_dna_id)]

#
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
  mdata_st3k = read.delim("stanford3k_metadata2.txt",sep=" ")
})

# Same data from Kirstie
mdata_kirs = read.xlsx2("Stanford_Ashley_MEGAv2_n3484_DNAReport_Kirstie_dw.xlsx",1)
try({mdata_kirs = read.delim("metadata/Stanford_Ashley_MEGAv2_n3484_DNAReport_Kirstie_dw.tsv")})

################# New script mapping sample ids to their idat files ##################
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
# solve situations in which we have duplicated records
remove_dup_results_in_mapping_vector<-function(v,which_to_take="last"){
  newv = c()
  tt = table(v)
  uniques = names(which(tt==1))
  is_unique = is.element(v,set=uniques)
  newv[names(v)[is_unique]] = v[is_unique]
  non_u = names(which(tt>1))
  for(id in non_u){
    inds = which(v==id)
    if(which_to_take=="last"){
      inds = inds[length(inds)]
    }
    else{
      inds = inds[1]
    }
    newv[names(v)[inds]] = v[inds]
  }
  return(newv)
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
alternate_id2sample_id = reverse_names(sample_id_to_alternate_id)
alternate_id2sample_id = c(alternate_id2sample_id,clinical_id2sample_id)
length(alternate_id2sample_id)

# use the stanford3k file to map dna ids to sample ids
dnaid2sampleid_st3k = as.character(mdata_st3k$FID)
names(dnaid2sampleid_st3k) = mdata_st3k$IID
# use Kirsties' metadata file to map dna ids to sample ids
dnaid2sampleid_kir = as.character(mdata_kirs$id2)
names(dnaid2sampleid_kir) = mdata_kirs$DNA_ID
length(dnaid2sampleid_kir)
length(dnaid2sampleid_st3k)
length(intersect(dnaid2sampleid_st3k,dnaid2sampleid_st3k)) # roughly the same

# This is the main table we create here. The addition rules below are ordered.
sample_id2idat_file = c()
unmapped_observed_idats = c()

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
print(dim(sample_id2idat_file))
length(unique(rownames(sample_id2idat_file)))

# Rule 1.1: use Mikael's CL id mappings from June 2018
new_cl_info = read.xlsx2('fitness_genetics_ashleylab_metadata.xlsx',"june2018_mikael_cl_id_mapping",stringsAsFactors = F)
correct_elite_ids<-function(x){
  arr = strsplit(x,split='_')[[1]]
  if(length(arr)<=1){return(x)}
  arr = arr[1:2]
  if(nchar(arr[1])==1){
    arr[1] = paste("0",arr[1],sep="")
  }
  newx = paste(arr,collapse="_")
}
new_cl_mapping = new_cl_info[,2]
names(new_cl_mapping) = sapply(new_cl_info[,1],correct_elite_ids)
new_cl_mapping = reverse_names(new_cl_mapping)
cl_ids_from_dna_ids = unlist(sapply(idats_dna_id,function(x){arr=strsplit(x,split='_')[[1]];arr[length(arr)]}))
new_cls_with_dnaid = intersect(names(new_cl_mapping),cl_ids_from_dna_ids)
new_cls_with_dnaid_sample_ids = new_cl_mapping[new_cls_with_dnaid]
new_cls_with_dnaid_sample_ids = new_cls_with_dnaid_sample_ids[
  !is.element(new_cls_with_dnaid_sample_ids,rownames(sample_id2idat_file))]
length(new_cls_with_dnaid_sample_ids)
for(cl in names(new_cls_with_dnaid_sample_ids)){
  s = new_cls_with_dnaid_sample_ids[cl]
  ind = which(cl_ids_from_dna_ids==cl)
  curr_barcodes = idat_barcodes[ind]
  curr_locs = idat_locs[ind]
  sample_id2idat_file = rbind(sample_id2idat_file,c(curr_barcodes[1],curr_locs[1]))
  rownames(sample_id2idat_file)[nrow(sample_id2idat_file)] = s
}
print(dim(sample_id2idat_file))
length(unique(rownames(sample_id2idat_file)))

# Rule 2: use the DNA IDs in the mainfile to map to the idat files
dnaid2sampleid_mainfile = reverse_names(sample_id2dna_id_mainfile)
mainfile_shared_dnaids_with_idats = intersect(names(dnaid2sampleid_mainfile),names(idat_barcodes))
curr_dnaid2sample_id = dnaid2sampleid_mainfile[mainfile_shared_dnaids_with_idats]
print(length(curr_dnaid2sample_id))
print(length(unique(curr_dnaid2sample_id)))
curr_dnaid2sample_id = remove_dup_results_in_mapping_vector(curr_dnaid2sample_id)
print(length(curr_dnaid2sample_id))
print(length(unique(curr_dnaid2sample_id)))

# As of may 2018 this is empty:
unmapped_but_in_mainfile = names(which(""==curr_dnaid2sample_id | is.na(curr_dnaid2sample_id)))
print(dim(sample_id2idat_file))
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,curr_dnaid2sample_id,idat_barcodes,idat_locs)
print(dim(sample_id2idat_file))
length(unique(rownames(sample_id2idat_file)))
# table(is.element(rownames(sample_id2idat_file),set=all_sample_data$Sample_ID))
# table(is.element(rownames(sample_id2idat_file),set=names(alternate_id2sample_id)))

# COMMENT: as of June 12 2018: Rules 3 and 4 do not add samples

# Rule 3: use stanford 3k DNA IDs to map to the idat files
st3k_shared_dnaids_with_idats = intersect(names(dnaid2sampleid_st3k),names(idat_barcodes))
curr_dnaid2sample_id = dnaid2sampleid_st3k[st3k_shared_dnaids_with_idats]
table(is.element(curr_dnaid2sample_id,set=names(alternate_id2sample_id)))
table(is.element(curr_dnaid2sample_id,set=all_sample_data$Sample_ID))
curr_dnaid2sample_id = curr_dnaid2sample_id[is.element(curr_dnaid2sample_id,set=all_sample_data$Sample_ID)]
print(dim(sample_id2idat_file))
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,curr_dnaid2sample_id,idat_barcodes,idat_locs)
print(dim(sample_id2idat_file))
table(is.element(rownames(sample_id2idat_file),set=all_sample_data$Sample_ID))
table(is.element(rownames(sample_id2idat_file),set=names(sample_id_to_alternate_id)))

# Rule 4: use the dna ids in the idat file names
sample_ids_from_dna_ids = unlist(sapply(idats_dna_id,function(x){arr=strsplit(x,split='_')[[1]];arr[length(arr)]}))
sample_ids_from_dna_ids = sample_ids_from_dna_ids[is.element(sample_ids_from_dna_ids,
                                                             set=union(all_sample_data$Sample_ID,names(alternate_id2sample_id)))]
curr_alt_ids = is.element(sample_ids_from_dna_ids,set=names(alternate_id2sample_id))
table(curr_alt_ids)
print(dim(sample_id2idat_file))
sample_id2idat_file = add_based_on_new_ids(sample_id2idat_file,sample_ids_from_dna_ids,idat_barcodes,idat_locs)
print(dim(sample_id2idat_file))

# add phenotypes and other features
m = c()
for (s in rownames(sample_id2idat_file)){
  ind = which(all_sample_data$Sample_ID==s)[1]
  m = rbind(m,all_sample_data[ind,])
}
rownames(m) = rownames(sample_id2idat_file)
sample_id2idat_file = cbind(sample_id2idat_file,m)
table(sample_id2idat_file$Cohort)

# print the file
# Version 1:
write.table(sample_id2idat_file,file="merged_metadata_file_stanford3k_elite_cooper.txt",quote=F,sep="\t")
# Jan 2019: version 2 after correcting some mis-formatted sample ids
write.table(sample_id2idat_file,file="merged_metadata_file_stanford3k_elite_cooper_v2_jan_2019.txt",quote=F,sep="\t")
is.element("62E02",rownames(sample_id2idat_file))

# # QC
# # compare the two file versions: befor and after Jan 2019:
# x1 = read.delim("merged_metadata_file_stanford3k_elite_cooper.txt",row.names = 1,stringsAsFactors = F)
# x2 = read.delim("merged_metadata_file_stanford3k_elite_cooper_v2_jan_2019.txt",row.names = 1,stringsAsFactors = F)
# setdiff(rownames(x2),rownames(x1))
# inds = intersect(rownames(x2),rownames(x1))
# all(x1[inds,]==x2[inds,],na.rm=T)
# all(which(is.na(x1))==which(is.na(x2)))
# sort(colSums(x1!=x2,na.rm=T))

# # Get stats by cohorts
# curr_rows = is.element(all_sample_data$Sample_ID,set=rownames(sample_id2idat_file))
# table(all_sample_data$Cohort[curr_rows])
# table(table(sample_id2idat_file[,1]))
# length(setdiff(unique(idat_barcodes),unique(sample_id2idat_file[,1])))
# length(intersect(unique(idat_barcodes),unique(sample_id2idat_file[,1])))
# table(is.element(rownames(sample_id2idat_file),set=sample_id_to_alternate_id))

# # create a sample sheet for illumina
# mdata_batch1 = as.matrix(read.csv('1092_samples_report_1.csv',skip = 15,header = T))
# mdata_batch2 = as.matrix(read.csv('485_samples_report_2.csv',skip = 15,header = T))
# 
# m = sample_id2idat_file
# dna_ids1 = sample_id2dna_id_mainfile[rownames(m)]
# dna_ids1 = sample_id2dna_id[rownames(m)]
# 
# # Additional QC tests
# # Compare sample id to dna id mapping
# s3k_malene_intersect = intersect(mdata_elite$Sample_ID,mdata_st3k$FID)
# table(is.element(s3k_malene_intersect,set=rownames(mdata_reports)))
# x1 = mdata_elite[is.element(mdata_elite$Sample_ID,set=s3k_malene_intersect),]
# x2 = mdata_st3k[is.element(mdata_st3k$FID,set=s3k_malene_intersect),]
# discrepancies = c()
# for(i in 1:nrow(x1)){
#   id1 = as.character(x1$Sample_ID[i])
#   dnaid1 = as.character(x1$DNA_ID[i])
#   inds = which(x2$FID==id1)
#   dnaid2 = as.character(x2$IID[inds])
#   if(any(dnaid2!=dnaid1)){
#     print(id1)
#     discrepancies = rbind(discrepancies,c(id1,dnaid1,dnaid2[1]))
#   }
# }
# colnames(discrepancies) = c("SampleID","ELITE_sheet1","Stanford3k_sheet")
# write.csv(discrepancies,file="diff_elite_vs_stanford3k.csv")
# discrepancies[which(is.element(discrepancies[,1],set=rownames(mdata_reports))),]
