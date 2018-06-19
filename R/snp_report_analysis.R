
# This script goes over the snp report and provides some statistics
# Thresholds and parameters are based on:
# https://github.com/MerrimanLab/merrimanlab.github.io/wiki/Genome-Studio-Genotype-Quality-Control-Protocol

# in sherlock
snp_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
# locally
snp_report_file = "/Users/David/Desktop/elite/reports/no_reclustering_SNP_Table.txt"

# in sherlock
snps_to_remove_external = '/oak/stanford/groups/euan/projects/mega-cooper-elite-udn/data/snps_to_remove1.txt'
# locally
snps_to_remove_external = "/Users/David/Desktop/elite/reports/snps_to_remove1.txt"

# in sherlock 
sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
# locally
sample_report_file = "/Users/David/Desktop/elite/reports/no_reclustering_Samples_Table.txt"

# in sherlock 
processed_sample_report_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/integrated_sample_metadata_and_covariates.txt"
# locally
processed_sample_report_file = "/Users/David/Desktop/elite/analysis/integrated_sample_metadata_and_covariates.txt"

# in sherlock 
#old_recall_info = "/Users/David/Desktop/elite/metadata/1092_samples_report_table.tsv"
# locally
old_recall_info = "/Users/David/Desktop/elite/metadata/1092_samples_report_table.tsv"

# in sherlock 
#metadata_file = "/Users/David/Desktop/elite/metadata/1092_samples_report_table.tsv"
# locally
metadata_file = "/Users/David/Desktop/elite/metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt"

# in sherlock 
#gs_excluded_file = "/Users/David/Desktop/elite/metadata/1092_samples_report_table.tsv"
# locally
gs_excluded_file = "/Users/David/Desktop/elite/metadata/june_2018_integrated_info/low_call_rate_samples_genomestudio.txt"

# load SNP info
snp_data = read.delim(snp_report_file)
rownames(snp_data) = snp_data$Name
snps_to_remove_external_raw = read.delim(snps_to_remove_external,stringsAsFactors = F)
snps_to_remove_external = read.delim(snps_to_remove_external,stringsAsFactors = F)[,1]
length(snps_to_remove_external)
length(intersect(snp_data$Name,snps_to_remove_external))
setdiff(snps_to_remove_external,snp_data$Name)

# some simple tests on a random selection
snp_data_sample = snp_data[sample(1:nrow(snp_data))[1:100000],]
snp_data_excluded_snps = snp_data[intersect(snp_data$Name,snps_to_remove_external),]

quantile(snp_data_sample$ChiTest100)
quantile(snp_data_excluded_snps$ChiTest100)

quantile(snp_data_sample$Call.Freq)
quantile(snp_data_excluded_snps$Call.Freq)

quantile(snp_data_sample$Multi.EthnicGlobal_D1.bpm.Cluster.Sep)
quantile(snp_data_excluded_snps$Multi.EthnicGlobal_D1.bpm.Cluster.Sep)

quantile(snp_data_sample$Het.Excess)
quantile(snp_data_excluded_snps$Het.Excess)

# Useful stats
table(snp_data$Call.Freq<0.95)
table(snp_data$Call.Freq<0.95)/nrow(snp_data)
table(snp_data[snp_data$Call.Freq<0.93,]$Chr)
table(snp_data$Chr)

# # get a final set of excluded snps and compare
# snp_min_clustersep_thr = 0.5
# snp_min_call_rate = 0.95
# snp_min_het_ex = -0.3
# snp_max_het_ex = 0.2
# min_maf = 0.005
# snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
# snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
#   snp_data$Multi.EthnicGlobal_D1.bpm.Cluster.Sep < snp_min_clustersep_thr |
#   snp_data$Het.Excess < snp_min_het_ex |
#   snp_data$Het.Excess > snp_max_het_ex
# snps_to_exclude = snps_to_exclude & snp_data_autosomal_rows
# our_excluded_snps = snp_data$Name[snps_to_exclude]
# length(intersect(our_excluded_snps,snps_to_remove_external))
# table(snp_data$Minor.Freq < min_maf)

# Analysis of the sample data
sample_data = read.delim(sample_report_file,stringsAsFactors = F)
processed_sample_data = read.delim(processed_sample_report_file,stringsAsFactors = F)
old_recall_info = read.delim(old_recall_info)

old_call_rates1 = old_recall_info$RAW.DATA.CALL.RATE
old_call_rates2 = old_recall_info$CALL.RATE...CALL.FRAG...95....snps
names(old_call_rates1) = paste(old_recall_info$Serial.number.of.chips,old_recall_info$Sentri.x.Position,sep="_")
names(old_call_rates2) = names(old_call_rates1)

our_call_rates = sample_data$Call.Rate
names(our_call_rates) = paste(sample_data$Array.Info.Sentrix.ID,sample_data$Array.Info.Sentrix.Position,sep="_")

our_call_rates_proc = processed_sample_data$call_rates_after_filters
names(our_call_rates_proc) = rownames(processed_sample_data)

par(mfrow=c(1,2))
inds = intersect(names(old_call_rates1),names(our_call_rates))
plot(our_call_rates[inds],old_call_rates1[inds],pch=3,
     ylab = "Illu report call rate", xlab = "Our call rate",main="Raw data call rates");abline(0,1,lty=2,lwd=2,col="blue")

inds = intersect(names(old_call_rates2),names(our_call_rates_proc))
plot(our_call_rates_proc[inds],old_call_rates2[inds],pch=3,
     ylab = "Illu report call rate", xlab = "Our call rate",main="Filtered data call rates");abline(0,1,lty=2,lwd=2,col="blue")

plot(old_call_rates1,old_call_rates2,xlim=c(0.8,1),ylim=c(0.8,1));abline(0,1,lty=2,lwd=2,col="blue")

# Look at the samples excluded because they had low raw call rates in genomestudio
proj_metadata = read.delim(metadata_file)
gs_excluded = read.delim(gs_excluded_file)
write.table(
  cbind(gs_excluded$Sample.ID,
  proj_metadata[gs_excluded$Sample.ID,c("Sample_ID","Cohort")],
  gs_excluded$Call.Rate),
            file ="/Users/David/Desktop/elite/metadata/june_2018_integrated_info/low_call_rate_samples_genomestudio_our_ids.txt",
            sep="\t",quote=F)


