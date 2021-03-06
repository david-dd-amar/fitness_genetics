
# Number of controls for selecion:
N = 20000

# Input files: locally
# Metadata: an excel file with multiple sheets. One per project.
our_metadata_file = "/Users/David/Desktop/elite/metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt"
# Ukbb phenotypic data
ukbb_pheno_file_raw = "/Users/David/Desktop/ukbb/biobank_collated_pheno_data.RData"
ukbb_pheno_file_processed = "/Users/David/Desktop/ukbb/covariate_matrix.RData"
# Load euro samples
euro_sample_ids_file = "/Users/David/Desktop/ukbb/pca_results_v2_chrom1_euro.eigenvec"
# Load genotypes samples after qc
genotyped_samples = "/Users/David/Desktop/elite/ukbb/ukb_imp_chr1_v2.fam"
# Subjects that were found fit for the exercise test
ukbb_exercise_subjects = "/Users/David/Desktop/ukbb/gwas/simple_norm/RHR_fitness_subjects.txt"

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
eu_ids = read.delim(euro_sample_ids_file)
eu_ids = as.character(eu_ids[,1])
gc()
geno_samples = intersect(geno_samples,eu_ids)

our_metadata = read.delim(our_metadata_file,stringsAsFactors = F)
our_metadata = our_metadata[is.element(our_metadata$Cohort,set=c("ELITE","Cooper")),]
hist(as.numeric(our_metadata$Age..at.test.))
table(our_metadata$Sex_input_data)/nrow(our_metadata)

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
ukbb_exercise_subjects = read.table(ukbb_exercise_subjects,header = T)
ukbb_exercise_subjects = as.character(ukbb_exercise_subjects[,1])
hist(ages)
hist(pheno_data[ukbb_exercise_subjects,defined_features["age","ukbb"]])
# select subjects aged up to 60
geno_samples = geno_samples[ages <=60]
length(geno_samples)
geno_samples = intersect(geno_samples,ukbb_exercise_subjects)
length(geno_samples)

# Select based on sex only
our_sex_dist = table(our_metadata$Sex_input_data)
male_prob = (our_sex_dist / nrow(our_metadata))["M"]
geno_samples_sex = pheno_data[geno_samples,defined_features["sex","ukbb"]]
geno_samples_sex_male = geno_samples_sex=="Male"
N_male = N*male_prob
male_samples = geno_samples[sample(which(geno_samples_sex_male))[1:N_male]]
N_female = N - length(male_samples)
female_samples = geno_samples[sample(which(!geno_samples_sex_male))[1:N_female]]

selected_samples = sort(c(male_samples,female_samples))
length(selected_samples)
length(male_samples)
write.table(cbind(selected_samples,selected_samples),file="/Users/David/Desktop/elite/ukbb/tmp_20k_rand_controls_sex_age.txt",
            sep="\t",row.names = F,col.names = F,quote=F)
selected_samples = cbind(selected_samples,
                         pheno_data[selected_samples,defined_features["sex","ukbb"]],
                         pheno_data[selected_samples,defined_features["age","ukbb"]])
colnames(selected_samples) = c("ID","sex","age")
write.table(selected_samples,file="/Users/David/Desktop/elite/ukbb/tmp_20k_rand_controls_sex_age_with_info2.txt",
            sep="\t",row.names = F,col.names = F,quote=F)
