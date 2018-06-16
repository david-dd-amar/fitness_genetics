
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
snp_data = read.delim(snp_report_file)
rownames(snp_data) = snp_data$Name
snps_to_remove_external_raw = read.delim(snps_to_remove_external,stringsAsFactors = F)
snps_to_remove_external = read.delim(snps_to_remove_external,stringsAsFactors = F)[,1]
length(snps_to_remove_external)
length(intersect(snp_data$Name,snps_to_remove_external))
setdiff(snps_to_remove_external,snp_data$Name)

table(grepl("^X",snps_to_remove_external))
table(grepl("^rs",snps_to_remove_external))

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

table(snp_data$Call.Freq<0.95)
table(snp_data$Call.Freq<0.95)/nrow(snp_data)
table(snp_data[snp_data$Call.Freq<0.93,]$Chr)
table(snp_data$Chr)

# get a final set of excluded snps and compare
snp_min_clustersep_thr = 0.5
snp_min_call_rate = 0.95
snp_min_het_ex = -0.3
snp_max_het_ex = 0.2
min_maf = 0.005
snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
  snp_data$Multi.EthnicGlobal_D1.bpm.Cluster.Sep < snp_min_clustersep_thr |
  snp_data$Het.Excess < snp_min_het_ex |
  snp_data$Het.Excess > snp_max_het_ex
snps_to_exclude = snps_to_exclude & snp_data_autosomal_rows
our_excluded_snps = snp_data$Name[snps_to_exclude]
length(intersect(our_excluded_snps,snps_to_remove_external))
table(snp_data$Minor.Freq < min_maf)

#####
sample_data = read.delim(sample_report_file,stringsAsFactors = F)

table(sample_data$Call.Rate<0.945)
table(sample_data$Call.Rate<0.98)
quantile(sample_data$Call.Rate)

hist(snp_data$Call.Freq,breaks = 100)

## Compare our call rates to those from Illumina
old_recall_info = "/Users/David/Desktop/elite/metadata/1092_samples_report_table.tsv"
old_recall_info = read.delim(old_recall_info)
our_call_rates = sample_data$Call.Rate
names(our_call_rates) = paste(sample_data$Array.Info.Sentrix.ID,sample_data$Array.Info.Sentrix.Position,sep="_")
old_call_rates = old_recall_info$RAW.DATA.CALL.RATE
old_call_rates = old_recall_info$CALL.RATE...CALL.FRAG...95....snps
names(old_call_rates) = paste(old_recall_info$Serial.number.of.chips,old_recall_info$Sentri.x.Position,sep="_")
inds = intersect(names(old_call_rates),names(our_call_rates))
plot(our_call_rates[inds],old_call_rates[inds]);abline(0,1)
hist(our_call_rates[inds] - old_call_rates[inds])
wilcox.test(our_call_rates[inds],old_call_rates[inds])






