# Input files for our analyses below:
# 1. "merged_metadata_mega.tsv" - the main non-genetic phenotypic data sheet, has IIDs as first column
# 2. mega_all_pca.phe - PCA results when using all samples, includes and "is_eu" column
# 3. merged_mega_data.eu.eigenvec - PCA of the EU subjects
# 4. Sex check files (for imputation)
library(data.table)

# 1. Read the main non-genetic data, reformat the sheepment data, define batch, and binary variables for GWAS
d = read.delim("merged_metadata_mega.tsv",stringsAsFactors=F,header=T)
#head(d)
d[d==""] = NA
d[d==" "] = NA
d[d=="\t"] = NA
rownames(d) = d[[1]]
#print(names(which(table(d[[1]])>1)))

# take care of shipment date
d[is.na(d[["Shipment.date"]]),"Shipment.date"] = "Other"
d[["Shipment.date"]] = gsub(" ","",d[["Shipment.date"]])
d[["Shipment.date"]] = gsub("\\)","",d[["Shipment.date"]])
d[["Shipment.date"]] = gsub("\\(","",d[["Shipment.date"]])
# table(d[["Shipment.date"]])
d$batch = as.numeric(as.factor(d[["Shipment.date"]]))

# take care of cohort
print("After reading the non-genetic data, these are the cohort numbers:")
print(table(d[["Cohort"]]))
# add all pairwise comparisons
# Cooper    ELITE genepool
d$cooper_vs_elite = NA
d$cooper_vs_elite[d[["Cohort"]]=="ELITE"] = 1
d$cooper_vs_elite[d[["Cohort"]]=="Cooper"] = 2
print(table(d$cooper_vs_elite))
d$cooper_vs_gp = NA
d$cooper_vs_gp[d[["Cohort"]]=="Cooper"] = 2
d$cooper_vs_gp[d[["Cohort"]]=="genepool"] = 1
print(table(d$cooper_vs_gp))
d$elite_vs_gp = NA
d$elite_vs_gp[d[["Cohort"]]=="genepool"] = 1
d$elite_vs_gp[d[["Cohort"]]=="ELITE"] = 2
print(table(d$elite_vs_gp))

selected_cols = c(
"mega_plink_id","Sample_ID","Shipment.date","Cohort","Sex_input_data","Age..at.test.","VO2max..ml.kg.min.","VO2max..l.","Treadmill.time.sec",
"Ethnicity","Weight_kg","Height_cm","cooper_vs_elite","cooper_vs_gp","elite_vs_gp"
)
d = d[,selected_cols]
colnames(d)[1] = "#IID"

print("Before adding PCs, dim is:")
print(dim(d))

# 2. Add the PCs and the EU IDs
d2 = read.delim("mega_all_pca.phe",stringsAsFactors=F,header=T)
d2 = d2[names(d2)!="cohort"]
#head(d2)
iid_col_d2 = names(d2)[grepl("IID",names(d2))]
d = merge(d,d2,by.x="#IID",by.y=iid_col_d2,all.x=T)
print("After adding PCs, dim is:")
print(dim(d))
print("Num NAs in PC1:")
print(sum(is.na(d$PC1)))
rownames(d) = d[[1]]

# 3. Add the EU PCs
d3 = read.delim("merged_mega_data.eu.eigenvec",stringsAsFactors=F,header=T)
colnames(d3) = gsub("PC","EU_PC",colnames(d3))
iid_col_d3 = names(d3)[grepl("IID",names(d3))]
d = merge(d,d3,by.x="#IID",by.y=iid_col_d3,all.x=T)
#head(d)
print("After  adding EU PCs, dim is:")
print(dim(d))
print("Num NAs in EU_PC1:")
print(sum(is.na(d$EU_PC1)))
print("Variable names after merges:")
print(names(d))
rownames(d) = d[[1]]

#4. Take care of sex:
pred1 = fread("/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/global_b37_no_reclustering/PLINK_020520_0243/global_b37_no_reclustering.sexcheck",
	stringsAsFactors=F,data.table=F)
sex1 = pred1$SNPSEX
names(sex1) = pred1$IID
pred2 = fread("/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/genomestudio_projects/consortium_b37_with_reclustering/PLINK_100520_1057/consortium_b37_with_reclustering.sexcheck",
	stringsAsFactors=F,data.table=F)
sex2 = pred2$SNPSEX
names(sex2) = pred2$IID
geno_predicted_sex = c(sex1,sex2)
geno_predicted_sex[geno_predicted_sex==0] = NA
d$geno_predicted_sex = NA
intr = intersect(names(geno_predicted_sex),d[[1]])
d[intr,"geno_predicted_sex"] = geno_predicted_sex[intr]
print("Agreement table between non-genetic data sex records and the genetic prediction")
print(table(d$Sex_input_data,d$geno_predicted_sex))

# Create the final sex column for plink (1=M,2=F)
d$sex = d$Sex_input_data
d$sex[d$sex == "M"] = 1
d$sex[d$sex == "F"] = 2
d$sex[is.na(d$sex)] = d$geno_predicted_sex[is.na(d$sex)]
print("After merging the two sources for sex info, these are the final numbers (1=M, taking the non-genetic info by default):")
print(table(d$sex))
print("Number of sex NAs:")
print(paste("Non genetic data:",sum(is.na(d$Sex_input_data))))
print(paste("Genetic data:",sum(is.na(d$geno_predicted_sex))))
print(paste("Merged:",sum(is.na(d$sex))))

write.table(d,file = "master_phe_mega.phe",row.names=F,col.names=T,quote=F,sep="\t")
