# Analysis of the GENEATHLETE dir
setwd("/oak/stanford/groups/euan/projects/fitness_genetics/GENEATHLETE/")
all_files = paste(getwd(),'/',list.files(),sep='')
# Check _1 vs _2
map_1_files = sort(all_files[grepl("_1.map",all_files)])
map_2_files = sort(all_files[grepl("_2.map",all_files)])
for(j in 1:length(map_1_files)){
  f1 = map_1_files[j]
  f2 = map_2_files[j]
  if(gsub("_\\d.map","",f1)!=gsub("_\\d.map","",f2)){print(paste("ERROR in file number:",j,f1,f2))}
  d1 = read.table(f1,header = F)
  d2 = read.table(f2,header = F)
  print(length(setdiff(d1$V2,d2$V2)))
  print(c(nrow(d1),nrow(d2)))
}