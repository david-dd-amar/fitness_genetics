
# Number of controls for selecion:
N = 10000

# Input files: locally
# Metadata: an excel file with multiple sheets. One per project.
our_metadata_file = "/Users/David/Desktop/elite/metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt"
# Ukbb phenotypic data
ukbb_pheno_file_raw = "/Users/David/Desktop/ukbb/biobank_collated_pheno_data.RData"
ukbb_pheno_file_processed = "/Users/David/Desktop/ukbb/covariate_matrix.RData"
# Load euro samples
euro_sample_ids_file = "/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec"
# Load genotypes samples after qc
genotyped_samples = "/Users/David/Desktop/elite/ukb_imp_chr1_v2.fam"

# # Input files: sherlock
# # Metadata: an excel file with multiple sheets. One per project.
# our_metadata_file = "/Users/David/Desktop/elite/fitness_genetics_ashleylab_metadata.xlsx"
# # Ukbb phenotypic data
# ukbb_pheno_file_raw = "/oak/stanford/groups/euan/projects/ukbb/data/bulk_data/ukb7454.r"
# # Load euro samples
# euro_sample_ids_file = "/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec"

# Load input data
load(ukbb_pheno_file_raw)
geno_samples = read.table(genotyped_samples,stringsAsFactors = F)
geno_samples = geno_samples[,1]
geno_samples = as.character(geno_samples)
gc()
# eu_ids = read.delim(euro_sample_ids_file)
# eu_ids = eu_ids[,1]
# gc()

our_metadata = read.delim(our_metadata_file)
our_metadata = our_metadata[is.element(our_metadata$Cohort,set=c("ELITE","Cooper")),]

# Define a set of features for the analysis
colnames(pheno_data)[grepl("Age when",colnames(pheno_data))]
defined_features = rbind(
  c("Age when attended assessment centre.0.0","Age..at.test."),
  c("Sex.0.0","Sex_input_data")
)
colnames(defined_features) = c("ukbb","our")
rownames(defined_features) = c("age","sex")

# Look at ukbb ages:
ages = pheno_data[geno_samples,defined_features["age","ukbb"]]
hist(ages)
# select subjects aged up to 50
geno_samples = geno_samples[ages <=50]
length(geno_samples)

# Attempt 1: select based on sex only
our_sex_dist = table(our_metadata$Sex_input_data)
male_prob = (our_sex_dist / nrow(our_metadata))["M"]
geno_samples_sex = pheno_data[geno_samples,defined_features["sex","ukbb"]]
geno_samples_sex_male = geno_samples_sex=="Male"
N_male = N*male_prob
male_samples = geno_samples[sample(which(geno_samples_sex_male))[1:N_male]]
N_female = N - length(male_samples)
female_samples = geno_samples[sample(which(!geno_samples_sex_male))[1:N_female]]

selected_samples = sort(c(male_samples,female_samples))
selected_samples = cbind(selected_samples,
                         pheno_data[selected_samples,defined_features["sex","ukbb"]],
                         pheno_data[selected_samples,defined_features["age","ukbb"]])

# Later: we can choose based on other features as well
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








