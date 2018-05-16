# Metadata: an excel file with multiple sheets. One per project.
metadata_file = "/Users/David/Desktop/elite/fitness_genetics_ashleylab_metadata.xlsx"

library(xlsx)
metadata_sheets = list()
for(j in 1:2){
  metadata_sheets[[j]] = read.xlsx2(metadata_file,j,stringsAsFactors=F)
}

# Ukbb phenotypic data
ukbb_pheno_file = "/Users/David/Desktop/ukbb/covariate_matrix.RData"
ukbb_pheno = get(load(ukbb_pheno_file))
# Load euro samples
euro_sample_ids_file = "/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec"
eu_ids = read.delim(euro_sample_ids_file)
eu_ids = eu_ids[,1]
gc()

# Define a set of features for the analysis
colnames(covariate_matrix)[grepl("Sex",colnames(covariate_matrix))]
defined_features = rbind(
  c("Age when attended assessment centre","Age..at.test."),
  c("Sex","Sex_input_data")
)
colnames(defined_features) = c("ukbb","our")
rownames(defined_features) = c("age","sex")

# Take a cohort from our data and get their control group
cohort_name = "elite"
mdata_sheet = metadata_sheets[[1]]
mdata_sheet = mdata_sheet[,-which(colnames(mdata_sheet)=="Comment")]
mdata_sheet = mdata_sheet[,-which(colnames(mdata_sheet)=="DNA_ID")]
mdata_sheet = unique(mdata_sheet)
cohort_rows = tolower(mdata_sheet$Cohort)==tolower(cohort_name)
table(cohort_rows)
cohort_mdata = mdata_sheet[cohort_rows,defined_features[,2]]
rownames(cohort_mdata) = mdata_sheet$Sample_ID[cohort_rows]
ukbb_mdata = covariate_matrix[,defined_features[,1]]

transform_mdata_to_numeric<-function(m,sexind=2,male_code = "Male"){
  m[,sexind] = as.numeric(m[,sexind]==male_code)
  m = as.matrix(m)
  mode(m) = "numeric"
  return(m)
}
cohort_mdata = transform_mdata_to_numeric(cohort_mdata,2,"M")
ukbb_mdata = transform_mdata_to_numeric(ukbb_mdata,2,"Male")

# Select controls by sample size






