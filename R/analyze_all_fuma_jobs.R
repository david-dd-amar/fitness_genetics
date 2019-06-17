setwd("/oak/stanford/groups/euan/projects/fitness_genetics/fuma_jobs/")

# Get all the .zip files and create a directory for each one with the results
allfiles = list.files("./",full.names = F)
zipfiles = allfiles[grepl(".zip$",allfiles)]
fuma_dirs = c()
for (zip in zipfiles){
  zipname = gsub(".zip","",zip)
  zipname = gsub("FUMA_job","",zipname)
  # system(paste("mkdir",zipname))
  # system(paste("unzip -o",zip,"-d",zipname))
  fuma_dirs = c(fuma_dirs,zipname)
}

# Read the job table (details per job)
metadata = read.table("job_table.txt",header=T,stringsAsFactors = F,row.names = 1)

# Thresholds for selecting 
CADD = 12.37
RDB = "^1"
FuncExclude = "intronic|intergenic"
R2 = 0.8

# Go over the jobs, get the variant table per job
library(data.table,lib = "~/R/packages")
for (job in rownames(metadata)){
  print("#####################################")
  print(metadata[job,1:2])
  
  if(length(list.files(job))==0){
    print("Job has no info:")
    print(job)
    print("Skipping")
    next
  }
  
  effects = fread(metadata[job,"Path"],data.table=F)
  snps = fread(paste(job,"/","snps.txt",sep=""),data.table=F)
  table(snps$func)
  # Exclude SNPs without p-value and "low" LD
  to_rem = is.na(snps$gwasP) & snps$r2 < R2
  to_rem[is.na(to_rem)] = F
  snps = snps[!to_rem,]
  # exclude variants with low CADD or low RDB
  to_rem = snps$CADD < CADD & grepl(FuncExclude,snps$func) & !grepl(RDB,snps$RDB)
  to_rem[is.na(to_rem) & grepl(FuncExclude,snps$func)] = T
  to_rem[is.na(to_rem)] = F
  snps = snps[!to_rem,]
  
  if(any(colnames(effects)=="OR") && any(colnames(effects)=="SNP")){
    effects = effects[is.element(effects[,"SNP"],set=snps[,"rsID"]) |
                        is.element(effects[,"BP"],set=snps[,"pos"]) ,]
  }
  else{
    effects = effects[is.element(effects[,3],set=snps[,"rsID"]) |
                        is.element(effects[,"POS"],set=snps[,"pos"]) ,]
  }
  annovar1 = fread(paste(job,"/","annov.txt",sep=""),data.table=F)
  annovar2 = fread(paste(job,"/","annov.stats.txt",sep=""),data.table = F)
  loci = fread(paste(job,"/","GenomicRiskLoci.txt",sep=""),data.table = F)
  # Catalog: filter to see the main results
  gcatalog = fread(paste(job,"/","gwascatalog.txt",sep=""),data.table = F)
  gcatalog = gcatalog[,c(1:5,13,31,34:35)]
  gcatalog = unique(gcatalog)
  # ANNOVAR: we are interested in exonic info here
  annovar1 = annovar1[!grepl(FuncExclude,annovar1$annot),]
  
  # Print the results to files
  print(paste("Creating results dir for job:",job))
  print(paste("Num loci:",nrow(loci)))
  currdir = paste(job,"_processed/",sep="")
  system(paste("mkdir",currdir))
  print(paste("SNP annot size:",c(nrow(snps))))
  write.table(snps,file=paste(currdir,"snps.txt",sep=""),
              row.names = F,col.names = T,quote=F,sep="\t")
  print(paste("ANNOVAR1 annot size:",c(nrow(annovar1))))
  write.table(annovar1,file=paste(currdir,"annovar1.txt",sep=""),
              row.names = F,col.names = T,quote=F,sep="\t")
  print(paste("gcatalog annot size:",c(nrow(gcatalog))))
  write.table(gcatalog,file=paste(currdir,"gcatalog.txt",sep=""),
              row.names = F,col.names = T,quote=F,sep="\t")
  print(paste("effects annot size:",c(nrow(effects))))
  write.table(effects,file=paste(currdir,"effects.txt",sep=""),
              row.names = F,col.names = T,quote=F,sep="\t")
}

# Add gene info to the results dirs
for (job in rownames(metadata)){
  system(paste("cp",paste(job,"/","genes.txt",sep=""),
               paste(job,"_processed/",sep="")))
}




