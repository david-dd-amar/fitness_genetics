library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
        print("usage:<cohort metadata file><reference metadata file><joint pca results file><output file name>")
	print("example for elite metadata: /oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt")
	print("example for ukbb metadata: /oak/stanford/groups/euan/projects/elite/ukbb_exome/pheno/ukb_pheno_slim.phe")
	print("example for pca results file: /oak/stanford/groups/euan/projects/elite/ukbb_exome/eu_beds/ukbb_elite_merged.filtered.eigenvec")
	print("assumption: in both the metadata files and in the pca results file, the second column is the unique sample id")
	print("the metadata files are mapped to their PCs, individuals without PCA results are excluded")
	print("this script is specific to our elite data analysis and assumes that the metadata files have the following columns:")
	print("   Sex_input_data")
	print("   Age..at.test.")
	print("   BMI")
	print("   Ethnicity")
	print("   Weight")
	print("   Height")
	print("the reference metadata is assumed to be from ukbb and has the Rivaslab fields:")
	print("   INI21001 - BMI")
	print("   INI50 - height")
	print("   INI21002 - weight")
	print("The script also does the following")
	print("   code missing values (NAs) using -9")
	print("   sex is transformed to 1 for Female and 2 for Male")
	print("   a class column is added with value 1 for UKBB and 2 for ELITE")
        q(save="no")
}

#out="/oak/stanford/groups/euan/projects/elite/ukbb_exome/pheno/"
#file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
## to save time we take a subset pf the master pheno file
#file2 = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phewas/resources/master.phe"
#d2 = fread(file2,stringsAsFactors=F,data.table=F)
##head(d2)
## INI21001 - BMI
## INI50 - height
## INI21002 - weight
#d2 = d2[,c(colnames(d2)[1:10],"INI21001","INI50","INI21002")]
#file2 = paste0(out,"ukb_pheno_slim.phe")
#fwrite(d2,file2,row.names=F,sep = " ",quote=F)

out="/oak/stanford/groups/euan/projects/elite/ukbb_exome/pheno/"
file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
file2 =	paste0(out,"ukb_pheno_slim.phe")
pcafile = "/oak/stanford/groups/euan/projects/elite/ukbb_exome/beds/ukbb_elite_merged.filtered.eigenvec"

out="/oak/stanford/groups/euan/projects/elite/ukbb_exome/eu_beds/"
file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
file2 = paste0("/oak/stanford/groups/euan/projects/elite/ukbb_exome/pheno/","ukb_pheno_slim.phe")
pcafile = "/oak/stanford/groups/euan/projects/elite/ukbb_exome/eu_beds/ukbb_elite_merged.filtered.eigenvec"

out = args[4]
file1 = args[1]
file2 = args[2]
pcafile = args[3]

d1 = fread(file1,stringsAsFactors=F,data.table=F)
#head(d1)
d2 = fread(file2,stringsAsFactors=F,data.table=F)
d2[d2==-9]=NA
#head(d2)
d3 = fread(pcafile,stringsAsFactors=F,data.table=F)
#head(d3)

#print(table(d3[,2] %in% d2[,2]))
#print(table(d3[,2] %in% d1[["Sample_ID"]]))

d1 = d1[d1[["Sample_ID"]] %in%  d3[,2],]
rownames(d1) = d1[["Sample_ID"]]
d2 = d2[d2[,2] %in% d3[,2],]
rownames(d2) = as.character(d2[,2])

merged_pheno = d3
merged_pheno$sex=NA
merged_pheno$age=NA
merged_pheno$class="UKB"
merged_pheno$VO2=NA
merged_pheno$height=NA
merged_pheno$weight=NA
merged_pheno$BMI=NA
merged_pheno$ethnicity=NA

rownames(merged_pheno) = as.character(merged_pheno[,2])
merged_pheno$class[is.na(as.numeric(merged_pheno[,2]))] = "ELITE"

# add BMI, VO2max..ml.kg.min., Ethnicity, Weight, Height
merged_pheno[rownames(d1),"sex"] = d1[,"Sex_input_data"] 
merged_pheno[rownames(d1),"age"] = d1[,"Age..at.test."]
merged_pheno[rownames(d1),"BMI"] = d1[,"BMI"]
merged_pheno[rownames(d1),"ethnicity"] = d1[,"Ethnicity"]
merged_pheno[rownames(d1),"weight"] = d1[,"Weight"]
merged_pheno[rownames(d1),"height"] = d1[,"Height"]
merged_pheno[rownames(d1),"VO2"] = d1[,"VO2max..ml.kg.min."]
print(table(merged_pheno[["sex"]]))

merged_pheno[rownames(d2),"sex"] = d2[,"sex"]
merged_pheno[rownames(d2),"age"] = d2[,"age"]
## INI21001 - BMI
## INI50 - height
## INI21002 - weight
merged_pheno[rownames(d2),"BMI"] = d2[,"INI21001"]
merged_pheno[rownames(d2),"height"] = d2[,"INI50"]
merged_pheno[rownames(d2),"weight"] = d2[,"INI21002"]
print(table(merged_pheno[["sex"]]))
print(table(merged_pheno[["class"]]))

# in ukbb - 0 is female 1 is male
merged_pheno[which(merged_pheno[,"sex"] == "0"),"sex"] = "F"
merged_pheno[which(merged_pheno[,"sex"] == "1"),"sex"] =	"M"
merged_pheno[which(merged_pheno[,"sex"] == "-9"),"sex"] =	NA
print(table(merged_pheno[["class"]],merged_pheno[["sex"]]))

print("looking at the merged dataset:")
#head(merged_pheno)

# convert to numeric
merged_pheno[["BMI"]] = as.numeric(as.character(merged_pheno[["BMI"]]))
merged_pheno[["weight"]] =	as.numeric(as.character(merged_pheno[["weight"]]))
merged_pheno[["height"]] =	as.numeric(as.character(merged_pheno[["height"]]))

sapply(merged_pheno,class)
table(is.na(merged_pheno))
merged_pheno[["weight"]][is.nan(merged_pheno[["weight"]])] = NA
merged_pheno[["height"]][is.nan(merged_pheno[["height"]])] = NA

quantile(merged_pheno[["weight"]],na.rm=T)
quantile(merged_pheno[["height"]],na.rm=T)

# TODO: make this work
#summary(lm(weight~factor(class),data=merged_pheno))
#summary(lm(BMI~factor(class),data=merged_pheno))
#summary(lm(height~factor(class),data=merged_pheno))

print("writing pheno file")
for(nn in names(merged_pheno)){
	merged_pheno[[nn]] = as.character(merged_pheno[[nn]])
}
print(table(merged_pheno==""))
merged_pheno[merged_pheno==""] = "-9"
merged_pheno[is.na(merged_pheno)] = "-9"
merged_pheno[merged_pheno=="NA"] = "-9"
# change binary variables to 0 and 1
merged_pheno[merged_pheno[["sex"]]=="F","sex"] = "1"
merged_pheno[merged_pheno[["sex"]]=="M","sex"] = "2"
merged_pheno[merged_pheno[["class"]]=="UKB","class"] = "1"
merged_pheno[merged_pheno[["class"]]=="ELITE","class"] = "2"
print(table(merged_pheno[["sex"]]))
print(table(merged_pheno[["class"]]))
fwrite(merged_pheno,file=out,sep="\t",quote=F,row.names=F,na="NA")
