dataset = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/impute2_1000gRef_out/"
chrs = paste("chr",c(1:22),sep="")
out_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool_imp/gwas_res/"
system(paste("mkdir",out_path))

# our annotation files
mega_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper_v2_jan_2019.txt"
genepool_covars_path = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/genepool/genepool_merged_info.txt"
sleep_study_annot_file = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/genepool/sleep_study_annots.csv"

out_anno_file = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper_with_gp.txt"

############################################################################
############################################################################
# Merge our metadata into one file
mega_anno = read.delim(mega_covars_path,stringsAsFactors = F)
gp_anno = read.delim(genepool_covars_path,stringsAsFactors = F,sep="\t")
gp_anno =  gp_anno["#N/A" != gp_anno[,2],]
# gp ages
gp_birth = as.numeric(sapply(gp_anno[,"BIRTH_DATE_SHIFTED"],function(x)strsplit(x,split="-")[[1]][3])) + 1900
gp_end_date = as.numeric(sapply(gp_anno[,"DEATH_DATE_SHIFTED"],function(x)strsplit(x,split="-")[[1]][3])) + 2000
gp_end_date[is.na(gp_end_date)] = 2016
gp_age = gp_end_date - gp_birth
gp_anno = cbind(gp_anno,gp_age)
gp_anno[gp_anno=="Missing"] = NA
# map ids
gp_anno = gp_anno[is.element(gp_anno[,2],set=rownames(mega_anno)),]
# solve duplications
dups = names(which(table(gp_anno[,2])>1))
to_rem = rep(F,nrow(gp_anno))
for(id in dups){
  curr_inds = which(gp_anno[,2]==id)
  num_nas = rowSums(is.na(gp_anno[curr_inds,]))
  curr_ind = curr_inds[num_nas==min(num_nas)][1]
  curr_inds = setdiff(curr_inds,curr_ind)
  to_rem[curr_inds]=T
}
gp_anno  = gp_anno[!to_rem,]
rownames(gp_anno) = gp_anno[,2]
# compare to our metadata
table(mega_anno[rownames(gp_anno),"Cohort"])
gp_anno[mega_anno[rownames(gp_anno),"Cohort"]!="genepool",]
table(mega_anno[rownames(gp_anno),"Sex_input_data" ],gp_anno[,"GENDER"])
table(mega_anno[rownames(gp_anno),"Sex_genotyped" ],gp_anno[,"GENDER"])

# Update our files:
mega_anno[rownames(gp_anno),"Sex_input_data" ] = gp_anno[,"GENDER"]
mega_anno[rownames(gp_anno),"Age..at.test." ] = gp_anno[,"gp_age"]

# Update the ethnicity col
colnames(mega_anno)[grepl("Ethnic",colnames(mega_anno))] = "Ethnicity"
mega_anno[rownames(gp_anno),"Ethnicity" ] = paste(gp_anno[,"ETHNICITY"],gp_anno[,"PRIMARY_RACE"],sep=";")

############################################################################
############################################################################
# Add in the sleep study information
sleep_anno = read.csv(sleep_study_annot_file,row.names = 1,stringsAsFactors = F)[,1:4]
# format the data nicely
sleep_anno[sleep_anno$Gender=="m","Gender"] = "Male"
sleep_anno[sleep_anno$Gender=="M","Gender"] = "Male"
sleep_anno[sleep_anno$Gender=="F","Gender"] = "Female"
sleep_anno[sleep_anno$Gender=="f","Gender"] = "Female"
sleep_anno[sleep_anno$Ethnicity=="Caucasian","Ethnicity"] = "Non-Hispanic;White"
# add the data into ours
inds = intersect(rownames(sleep_anno),rownames(mega_anno))
mega_anno[inds,"Sex_input_data" ] = sleep_anno[inds,"Gender"]
mega_anno[inds,"Age..at.test." ] = sleep_anno[inds,"Age"]
mega_anno[inds,"Ethnicity" ] = sleep_anno[inds,"Ethnicity"]

write.table(mega_anno,file=out_anno_file,row.names = T,col.names = T,sep="\t",quote=F)



