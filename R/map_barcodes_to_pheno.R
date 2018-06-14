library(xlsx)
fitness_metadata_sheets = "/Users/David/Desktop/elite/metadata/fitness_genetics_ashleylab_metadata.xlsx"
sample2codes_path = "/Users/David/Desktop/elite/metadata/merged_metadata_file_stanford3k_elite_cooper.csv"

all_metadata = read.xlsx(fitness_metadata_sheets,1)
failed_samples = !(is.na(all_metadata$Comment) | all_metadata$Comment!="Failed Sample")
failed_samples_data = all_metadata[failed_samples,]
#all_metadata = all_metadata[!failed_samples,]
all_codes = read.csv(sample2codes_path,stringsAsFactors = F)
table(all_metadata$Cohort)

plink_ids = paste(all_codes$SentrixBarcode_A,all_codes$SentrixPosition_A,sep="_")
our_ids = all_codes$Sample_ID

correct_numeric_id1<-function(x){
  num_zeros = 0
  j = nchar(x)
  while (j>1 && substr(x,j,j)=="0") {
    j = j-1
    num_zeros = num_zeros+1
  }
  if(num_zeros==0){return(x)}
  num_zeros = as.character(num_zeros)
  if(nchar(num_zeros)==1){
    num_zeros = paste("0",num_zeros,sep='')
  }
  newx = paste(
    substr(x,1,j),"E",num_zeros,sep=""
  )
  return(newx)
}
# correct_numeric_id("62")
# correct_numeric_id("620")
# correct_numeric_id("6200")
# correct_numeric_id("62000")
# correct_numeric_id("6200000000000")

get_id_row_from_spsheet<-function(id,all_metadata){
  curr_row=NA
  for(j in 1:3){
    curr_ind = which(all_metadata[,j]==id)[1]
    if(length(curr_ind)>0){
      curr_row = curr_ind
      break
    }
  }
  return(curr_row)
}

our_ids2metadata_row = c()
for(id in as.character(our_ids)){
  assign("last.warning", NULL, envir = baseenv())
  curr_row=get_id_row_from_spsheet(id,all_metadata)
  if(is.na(curr_row)){
    curr_row = get_id_row_from_spsheet(correct_numeric_id(id),all_metadata)
  }
  if(is.na(curr_row)){
    curr_row = get_id_row_from_spsheet(as.character(as.numeric(id)),all_metadata)
  }
  our_ids2metadata_row[as.character(id)] = curr_row
  if(!is.null(warnings())){
    print(id)
    break
  }
}
table(is.na(our_ids2metadata_row))
which(is.na(our_ids2metadata_row))
our_ids2metadata_row = our_ids2metadata_row[!is.na(our_ids2metadata_row)]

ourid2plink = plink_ids
names(ourid2plink) = our_ids

mapped_table = table(ourid2plink[names(our_ids2metadata_row)])
dup_codes = names(which(mapped_table>1))
our_ids[is.element(plink_ids,set=dup_codes)]

pheno_data_for_analysis = all_metadata[our_ids2metadata_row,]
rownames(pheno_data_for_analysis) = ourid2plink[names(our_ids2metadata_row)]

