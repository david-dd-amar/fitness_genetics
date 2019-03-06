# This script handles our entire GWAS analysis flow
# It creates a directory with all input, sh, log, err, and output files

####################################################################################################
####################################################################################################
####################################################################################################

# Define analysis parameters for the Illumina reports
autosomal_chrs = T
snp_min_clustersep_thr = 0.3
snp_min_call_rate = 0.95
snp_min_het_ex = -0.4
snp_max_het_ex = 0.4
min_maf = 0.001
run_loacally = F
num_pca_clusters = 3
# analysis_cohorts = "elite" # If null then use all cohorts
# analysis_cohorts = NULL # If null then use all cohorts

# Define analysis parameters for the Illumina reports
initial_subj_min_call_rate = 0.95
final_subj_min_call_rate = 0.98
snp_het_p = 1e-4

script_file = "/home/users/davidama/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

# Define Input parameters
# Original PLINK files
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/"
# ped_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.ped"
# map_file = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.map"
# New fwd strand files (August 2018)
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/"
# Alternative for elite only: August 2018
# job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/elite_only/our_prepro/"
# September 2018
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/"
# Dec 2018: add Genepool
job_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/"
# set the job's directory
try({system(paste("mkdir",job_dir),wait = T)})
setwd(job_dir)

# Each recalling has a set of parameters
# # 1. MEGG: old run with MEGC samples
# ped_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.ped"
# map_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering.map"
# input_bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_fwd_strand/recall_may_2018_without_reclustering"
# snp_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_SNP_Table.txt"
# sample_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/reports/no_reclustering_Samples_Table.txt"
# 1. MEGG: recalled alone (MEGG: MEGA global)
ped_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/raw.ped"
map_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/raw.map"
input_bfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/raw"
snp_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/SNP_Table.txt"
sample_report_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/Samples_Table.txt"
# 2. MEGC (MEGA Consortium)
ped_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data.ped"
map_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data.map"
input_bfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/plink_test_mega_consortium_data"
snp_report_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/SNP_Table.txt"
sample_report_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/Samples_Table.txt"
bad_snps_file = "/oak/stanford/groups/euan/projects/fitness_genetics/bad_mega_snps.txt"

# Our metadata
# Assumes that the first two columns are barcode and position. Third column is our sample ID.
sample_metadata = "/oak/stanford/groups/euan/projects/fitness_genetics/metadata/merged_metadata_file_stanford3k_elite_cooper.txt"
sample_metadata_raw = read.delim(sample_metadata,stringsAsFactors = F)
sample_metadata_raw = correct_dups_in_sample_metadata(sample_metadata_raw)

# # September 2018 we do not have genepool's metadata: we ignore this until we get it
# # This filter was later ignored during the Dec 2018 preprocessing
# sample_metadata_raw = sample_metadata_raw[sample_metadata_raw$Cohort!="genepool",]
# sample_metadata_raw = sample_metadata_raw[!is.na(sample_metadata_raw[,1]),]
# rownames(sample_metadata_raw) = apply(sample_metadata_raw[,1:2],1,paste,collapse="_")

# Dec 2018: get some stats for genepool samples
# gp_samples = sample_metadata_raw[sample_metadata_raw$Cohort=="genepool",]

####################################################################################################
####################################################################################################
####################################################################################################
# We need to merge the two MEGA sub datasets that were called separately
# Read the reports and compare
snp_data1 = read.delim(snp_report_file1)
sample_data1 = read.delim(sample_report_file1,stringsAsFactors = F)
rownames(sample_data1) = sample_data1$Sample.ID

snp_data2 = read.delim(snp_report_file2)
sample_data2 = read.delim(sample_report_file2,stringsAsFactors = F)
rownames(sample_data2) = sample_data2$Sample.ID

# Compare some stats
rownames(snp_data1) = snp_data1$Name
rownames(snp_data2) = snp_data2$Name
inds = intersect(rownames(snp_data1),rownames(snp_data2))
snp_samp = sample(inds)[1:2000]
cor(snp_data1[snp_samp,]$Call.Freq,snp_data2[snp_samp,]$Call.Freq,method="spearman")
cor(snp_data1[snp_samp,]$Het.Excess,snp_data2[snp_samp,]$Het.Excess,method="spearman")
snp_samp = snp_samp[!is.na(snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep)]
cor(snp_data1[snp_samp,]$Multi.EthnicGlobal_D1.bpm.Cluster.Sep,snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep,method="spearman")
table(snp_data2[snp_samp,]$MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep>0.3,snp_data1[snp_samp,]$Multi.EthnicGlobal_D1.bpm.Cluster.Sep>0.3)

# snp_data1["rs1747677",]

# some preprocessing of metadata
analyze_snp_report_get_snps_to_exclude<-function(snp_data,autosomal_chrs,snp_min_call_rate,
                                                 snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
                                                 clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep"){
  snp_data_autosomal_rows = grepl("^\\d+$",snp_data$Chr)
  snps_to_exclude =  snp_data$Call.Freq < snp_min_call_rate |
    snp_data[[clust_sep_col_name]] < snp_min_clustersep_thr |
    snp_data$Het.Excess < snp_min_het_ex |
    snp_data$Het.Excess > snp_max_het_ex
  if(autosomal_chrs){
    snps_to_exclude = snps_to_exclude & snp_data_autosomal_rows
  }
  if(!autosomal_chrs){
    snps_to_exclude = snps_to_exclude & !snp_data_autosomal_rows
  }
  return(snp_data$Name[snps_to_exclude])
}
excl1 = analyze_snp_report_get_snps_to_exclude(snp_data1,autosomal_chrs,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "Multi.EthnicGlobal_D1.bpm.Cluster.Sep")
excl2 = analyze_snp_report_get_snps_to_exclude(snp_data2,autosomal_chrs,snp_min_call_rate,
      snp_min_clustersep_thr,snp_min_het_ex,snp_max_het_ex,
      clust_sep_col_name = "MEGA_Consortium_v2_15070954_A2.bpm.Cluster.Sep")
excl3 = read.table(bad_snps_file,stringsAsFactors = F)[,1]

print("Initial QC using genomestudios reports")
print(paste("excluded from first file:",length(excl1)))
print(paste("excluded from secind file:",length(excl2)))
print(paste("intersect:",length(intersect(excl1,excl2))))

# Final snp set for subsequent analyses
snp_set = intersect(rownames(snp_data1),rownames(snp_data2))
snp_set = setdiff(snp_set,excl1)
snp_set = setdiff(snp_set,excl2)
snp_set = setdiff(snp_set,excl3)
write.table(t(t(snp_set)),paste(job_dir,"snp_set.txt",sep=""), row.names = F,col.names = F,quote = F)

# # Sanity check
# snp_set = readLines(paste(job_dir,"snp_set.txt",sep=""))
# bim1 = read.table(paste(input_bfile1,".bim",sep=""),stringsAsFactors = F)
# table(is.element(bim1[,2],set=snp_set))
# bim2 = read.table(paste(input_bfile1,".bim",sep=""),stringsAsFactors = F)
# table(is.element(bim2[,2],set=snp_set))

# take care of the FIDs
fam1 = read.table(paste(input_bfile1,".fam",sep=""),stringsAsFactors = F)
update_ids1 = cbind(fam1[,1:2],paste("file1_",fam1[,1],sep=""),fam1[,2])
write.table(file=paste(job_dir,"update_ids1.txt",sep=""),update_ids1,sep="\t",row.names = F,col.names = F,quote = F)
fam2 = read.table(paste(input_bfile2,".fam",sep=""),stringsAsFactors = F)
# for file 1 we want to remove samples that appear in file 2
inds = !is.element(fam1[,2],set=fam2[,2])
inds = inds & is.element(fam1[,2],set=rownames(sample_metadata_raw))
file1_samples_to_keep = fam1[inds,]
write.table(file=paste(job_dir,"file1_samples_to_keep.txt",sep=""),file1_samples_to_keep,sep="\t",row.names = F,col.names = F,quote = F)
write.table(file=paste(job_dir,"update_ids1.txt",sep=""),update_ids1[inds,],sep="\t",row.names = F,col.names = F,quote = F)

inds2 = is.element(fam2[,2],set=rownames(sample_metadata_raw))
file2_samples_to_keep = fam2[inds2,]
write.table(file=paste(job_dir,"file2_samples_to_keep.txt",sep=""),file2_samples_to_keep,sep="\t",row.names = F,col.names = F,quote = F)

# Correct File 2 before merge
# Remove samples that are not in the metadata
err_path = paste(job_dir,"prepare_bfile2.err",sep="")
log_path = paste(job_dir,"prepare_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile2,
                 "--keep",paste(job_dir,"file2_samples_to_keep.txt",sep=""),
                 "--extract",paste(job_dir,"snp_set.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "prepare_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))

# Correct File 1 before merge
# Remove samples that appear in file 2 (the consortium file)
err_path = paste(job_dir,"prepare_bfile1.err",sep="")
log_path = paste(job_dir,"prepare_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",input_bfile1,
                 "--keep",paste(job_dir,"file1_samples_to_keep.txt",sep=""),
                 "--extract",paste(job_dir,"snp_set.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "prepare_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

err_path = paste(job_dir,"update_ids_bfile1.err",sep="")
log_path = paste(job_dir,"update_ids_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--update-ids",paste(job_dir,"update_ids1.txt",sep=""),
                 "--make-bed --out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "update_ids_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job()

print("After initial qc, datas sizes are:")
print(paste("number of samples, file 1:",length(readLines(paste(job_dir,"bfile1.fam",sep="")))))
print(paste("number of snps, file 1:",length(readLines(paste(job_dir,"bfile1.bim",sep="")))))
print(paste("number of samples, file 2:",length(readLines(paste(job_dir,"bfile2.fam",sep="")))))
print(paste("number of snps, file 2:",length(readLines(paste(job_dir,"bfile2.bim",sep="")))))
# make sure there is no intersect
ids1 = read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F)[,2]
ids2 = read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of subjects in both files:",length(intersect(ids1,ids2))))

# Compare the bim files
bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F,row.names = 2)
bim2 = read.table(paste(job_dir,"bfile2.bim",sep=""),stringsAsFactors = F,row.names = 2)
all(dim(bim1)==dim(bim2))
length(intersect(rownames(bim1),rownames(bim2))) == nrow(bim1)
bim2 = bim2[rownames(bim1),]
table(bim1[,5]==bim2[,5])
same_snps = bim1[,5]==bim2[,5] & bim1[,4]==bim2[,4]
same_snps_diff_maf = bim1[,5]==bim2[,4] & bim1[,4]==bim2[,5]
zero_maf1 = bim1[,4]=="0"
zero_maf2 = bim2[,4]=="0"
maybe_rev_snps = !(same_snps | zero_maf2 | zero_maf1 | same_snps_diff_maf)
sum(maybe_rev_snps) == 0 # should be true

# Subject QC
# File 1
err_path = paste(job_dir,"subj_qc_bfile1.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--missing --het",
                 "--out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "subj_qc_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# File 2
err_path = paste(job_dir,"subj_qc_bfile2.err",sep="")
log_path = paste(job_dir,"subj_qc_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
                 "--missing --het",
                 "--out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "subj_qc_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# Analyze the results
# Subject call rates:
miss1 = read.table(paste(job_dir,"bfile1.imiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.imiss",sep=""),stringsAsFactors = F,header=T)
crs1 = 1-miss1$F_MISS;crs2 = 1-miss2$F_MISS
het1 = read.table(paste(job_dir,"bfile1.het",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.het",sep=""),stringsAsFactors = F,header=T)
# September 2018: we observed outlier subjects in file 2 only
to_rem2 = miss2[crs2<initial_subj_min_call_rate | het2$F< -0.25 | het2$F > 0.25 ,1:2]
to_rem1 = miss1[crs1<initial_subj_min_call_rate | het1$F< -0.25 | het1$F > 0.25 ,1:2]
remove_subjects_using_plink(paste(job_dir,"bfile2",sep=""),to_rem2,job_dir,"file2_initial_subj_qc","bfile2",
                                      batch_script_func=get_sh_default_prefix)
remove_subjects_using_plink(paste(job_dir,"bfile1",sep=""),to_rem1,job_dir,"file1_initial_subj_qc","bfile1",
                            batch_script_func=get_sh_default_prefix)
wait_for_job(60)
print("After initial qc, datas sizes are:")
print(paste("number of samples, file 1:",length(readLines(paste(job_dir,"bfile1.fam",sep="")))))
print(paste("number of snps, file 1:",length(readLines(paste(job_dir,"bfile1.bim",sep="")))))
ids1 = read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F)[,2]
ids2 = read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids1,"Cohort"]=="Cooper")))
print(paste("number of samples, file 2:",length(readLines(paste(job_dir,"bfile2.fam",sep="")))))
print(paste("number of snps, file 2:",length(readLines(paste(job_dir,"bfile2.bim",sep="")))))
print(paste("removed samples in previous step, by cohort",table(sample_metadata_raw[to_rem1[,2],"Cohort"])))
print(paste("removed samples in previous step, by cohort",table(sample_metadata_raw[to_rem2[,2],"Cohort"])))

# SNP QC: Get missingness, het, call rates etc for each file
# File 1
err_path = paste(job_dir,"snp_c_bfile1.err",sep="")
log_path = paste(job_dir,"snp_c_bfile1.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--freq 	--hardy --missing",
                 "--out",paste(job_dir,"bfile1",sep=''))
curr_sh_file = "snp_c_bfile1.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# File 2
err_path = paste(job_dir,"snp_c_bfile2.err",sep="")
log_path = paste(job_dir,"snp_c_bfile2.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
                 "--freq 	--hardy --missing",
                 "--out",paste(job_dir,"bfile2",sep=''))
curr_sh_file = "snp_c_bfile2.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

# Analyze the SNP QC results
miss1 = read.table(paste(job_dir,"bfile1.lmiss",sep=""),stringsAsFactors = F,header=T)
miss2 = read.table(paste(job_dir,"bfile2.lmiss",sep=""),stringsAsFactors = F,header=T)
frq1 = read.table(paste(job_dir,"bfile1.frq",sep=""),stringsAsFactors = F,header=T)
frq2 = read.table(paste(job_dir,"bfile2.frq",sep=""),stringsAsFactors = F,header=T)
het1 = read.table(paste(job_dir,"bfile1.hwe",sep=""),stringsAsFactors = F,header=T)
het2 = read.table(paste(job_dir,"bfile2.hwe",sep=""),stringsAsFactors = F,header=T)

to_rem_file1 = union(miss1$SNP[1-miss1$F_MISS<snp_min_call_rate],
                     frq1$SNP[frq1$MAF < min_maf])
to_rem_file1 = union(to_rem_file1,het1$SNP[het1$P < snp_het_p])
to_rem_file2 = union(miss2$SNP[1-miss2$F_MISS<snp_min_call_rate],
                     frq2$SNP[frq2$MAF < min_maf])
to_rem_file2 = union(to_rem_file2,het2$SNP[het2$P < snp_het_p])

# Some numbers
sum(1-miss1$F_MISS<snp_min_call_rate, na.rm=T)
sum(frq1$MAF < min_maf, na.rm = T)
sum(het1$P < snp_het_p, na.rm=T)
x1 = frq1$MAF;names(x1) = frq1$SNP
x2 = snp_data1$Minor.Freq;names(x2) = snp_data1$Name
x1[is.na(x1)] = 0
cor(x2[names(x1)],x1)
print(paste("Second QC using plink, removed from file 1:",length(to_rem_file1)))
print(paste("Second QC using plink, removed from file 2:",length(to_rem_file2)))
print(paste("Second QC, overlap between the two files:",length(intersect(to_rem_file2,to_rem_file1))))

to_rem_files_union = union(to_rem_file1,to_rem_file2)
bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
snps_to_keep = setdiff(bim1[,2],to_rem_files_union)
print(paste("remaining number of snps after second qc:",length(snps_to_keep)))

# Reduce the datasets before the merge
extract_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile1",
    batch_script_func=get_sh_default_prefix)
extract_snps_using_plink(paste(job_dir,"bfile2",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile2",
                         batch_script_func=get_sh_default_prefix)
wait_for_job()

print("After initial qc, datas sizes are:")
print(paste("number of samples, file 1:",length(readLines(paste(job_dir,"bfile1.fam",sep="")))))
print(paste("number of snps, file 1:",length(readLines(paste(job_dir,"bfile1.bim",sep="")))))
ids1 = read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F)[,2]
ids2 = read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids1,"Cohort"]=="Cooper")))
print(paste("number of samples, file 2:",length(readLines(paste(job_dir,"bfile2.fam",sep="")))))
print(paste("number of snps, file 2:",length(readLines(paste(job_dir,"bfile2.bim",sep="")))))
print(paste("removed samples in previous step, by cohort",table(sample_metadata_raw[to_rem1[,2],"Cohort"])))
print(paste("removed samples in previous step, by cohort",table(sample_metadata_raw[to_rem2[,2],"Cohort"])))

# Read the bim files: we may need to flip some snps
bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
bim2 = read.table(paste(job_dir,"bfile2.bim",sep=""),stringsAsFactors = F)
rownames(bim1) = bim1[,2];rownames(bim2) = bim2[,2]
all(bim1[,2]==bim2[,2])
bim1 = bim1[rownames(bim2),]

xx = cbind(bim1[,5:6],bim2[,5:6])
all_num_alleles = apply(xx,1,get_num_alleles)
one_allele_appears_as_two<-function(x){return(get_num_alleles(x)==2 & sum(x=="0")==2)}
is_one_allele_appears_as_two = apply(xx,1,one_allele_appears_as_two)
xx[is_one_allele_appears_as_two,]
snps_to_flip = rownames(bim1)[is_one_allele_appears_as_two | all_num_alleles>2]
print(paste("comparing the two bim files, these should be removed or flipped:",length(snps_to_flip)))
# Dec 2018: the snps discovered above have alleles marked as "0 0" in one of the files. Exclude and do not flip.
snps_to_keep = setdiff(snps_to_keep,snps_to_flip)
extract_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile1",
                         batch_script_func=get_sh_default_prefix)
extract_snps_using_plink(paste(job_dir,"bfile2",sep=""),snps_to_keep,job_dir,"final_snps_to_keep","bfile2",
                         batch_script_func=get_sh_default_prefix)
wait_for_job(60)

# repl = apply(xx[snps_to_flip,1:2],1,flip_snp_info)
# xx[snps_to_flip,1] = repl[1,]
# xx[snps_to_flip,2] = repl[2,]
# all_num_alleles2 = apply(xx,1,get_num_alleles)
# table(all_num_alleles2)
# all(is.na(bim1[all_num_alleles2==0,][1:10,])) # should be all NAs because there are no such snps
# flip_snps_using_plink(paste(job_dir,"bfile1",sep=""),snps_to_flip,job_dir,"final_snps_to_flip","bfile1",
#                          batch_script_func=get_sh_default_prefix)
# wait_for_job()

####################################################################################################
####################################################################################################
####################################################################################################
# # Run ld reports (useful for ldscore analysis)
# # File 1
# err_path = paste(job_dir,"ld_bfile1.err",sep="")
# log_path = paste(job_dir,"ld_bfile1.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
#                  "--r2 --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0.1 --freq",
#                  "--out",paste(job_dir,"bfile1",sep=''))
# curr_sh_file = "ld_bfile1.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# # File 2
# err_path = paste(job_dir,"ld_bfile2.err",sep="")
# log_path = paste(job_dir,"ld_bfile2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"bfile2",sep=''),
#                  "--r2 --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0.1 --freq",
#                  "--out",paste(job_dir,"bfile2",sep=''))
# curr_sh_file = "ld_bfile2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# 
# # Analyze the results
# bim1 = read.table(paste(job_dir,"bfile1.bim",sep=""),stringsAsFactors = F)
# bim2 = read.table(paste(job_dir,"bfile2.bim",sep=""),stringsAsFactors = F)
# rownames(bim1) = bim1[,2];rownames(bim2) = bim2[,2]
# 
# ld1 = read.table(paste(job_dir,"bfile1.ld",sep=""),stringsAsFactors = F,header=T)
# ld1_snps = table(c(ld1[,6],ld1[,3]))
# 
# res = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/with_ukbb_1000g_sanity/gwas/gwas_three_groups_linear.ExerciseGroup.glm.linear.adjusted",stringsAsFactors = F)
# res_top_snps = res[res[,3]<1e-50,2]
# ld1_snps[res_top_snps]

####################################################################################################
####################################################################################################
####################################################################################################
# Merge
err_path = paste(job_dir,"merge_plink.err",sep="")
log_path = paste(job_dir,"merge_plink.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"bfile1",sep=''),
                 "--bmerge",paste(job_dir,"bfile2",sep=''),
                 "--make-bed --out",paste(job_dir,"merged_mega_data",sep=''))
curr_sh_file = "merge_plink.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)
print("After merge, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data.bim",sep="")))))

# Check flipscan
# create a phe file
fam1 = as.matrix(read.table(paste(job_dir,"bfile1.fam",sep=""),stringsAsFactors = F))
fam2 = as.matrix(read.table(paste(job_dir,"bfile2.fam",sep=""),stringsAsFactors = F))
fam1 = cbind(fam1,rep("1",nrow(fam1)))
fam2 = cbind(fam2,rep("2",nrow(fam2)))
rownames(fam1)=NULL;colnames(fam1)=NULL
rownames(fam2)=NULL;colnames(fam2)=NULL
phe = rbind(fam1,fam2)
phe = phe[,c(1:2,7)]
colnames(phe) = c("FID","IID","mega")
write.table(phe,file=paste(job_dir,"flipscan_mega_type.txt",sep=""),
            sep=" ",quote=F,row.names = F,col.names = T)
# Run flipscan
err_path = paste(job_dir,"mega_flipscan.err",sep="")
log_path = paste(job_dir,"mega_flipscan.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 "--flip-scan --allow-no-sex",
                 "--pheno",paste(job_dir,"flipscan_mega_type.txt",sep=""),
                 "--pheno-name mega",
                 "--out",paste(job_dir,"merged_mega_data",sep=''))
curr_sh_file = "mega_flipscan.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(120)

flipscan_res = readLines(paste(job_dir,"merged_mega_data.flipscan",sep=""))
arrs = strsplit(flipscan_res[-1],split="\\s+")
names(arrs) = sapply(arrs,function(x)x[3])
table(sapply(arrs,length))
flipscan_failures = sapply(arrs,length) > 11
flipscan_failures = sapply(arrs[flipscan_failures],function(x)x[3])
# sapply(arrs[flipscan_failures],function(x)x[2])
bim = read.table(paste(job_dir,"merged_mega_data.bim",sep=""),stringsAsFactors = F)
snps_to_keep = setdiff(bim[,2],flipscan_failures)
print(paste("Flipscan check, number of variants to remove:",length(flipscan_failures)))
extract_snps_using_plink(paste(job_dir,"merged_mega_data",sep=""),snps_to_keep,job_dir,
                         "_final_snps_to_keep_after_flipscan","merged_mega_data",
                         batch_script_func=get_sh_default_prefix)
wait_for_job(60)

print("After flipscan, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data.bim",sep="")))))
ids = read.table(paste(job_dir,"merged_mega_data.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]!="Cooper")))
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Cooper")))

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# impute sex: each of the Mega datasets separately
run_plink_sex_check_x_chrom<-function(bfile,out_path,
                                      batch_script_func=get_sh_default_prefix,...){
  err_path = paste(job_dir,"impute_sex.err",sep="")
  log_path = paste(job_dir,"impute_sex.log",sep="")
  curr_cmd = paste("plink --bfile",bfile,
                   "--check-sex --chr 22-24 --out",bfile)
  curr_sh_file = paste(out_path,"impute_sex.sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file))
}
run_plink_sex_check_x_chrom(input_bfile1,job_dir)
run_plink_sex_check_x_chrom(input_bfile2,job_dir)
wait_for_job()

####################################################################################################
####################################################################################################
####################################################################################################
# keep the desired chrs (default: 0-22)
# also, exclude low qc score snps 
# in theory we can use the ped itself for inference and filtering
# however, here we assume we have the snp report that we can use
# above the list of excluded snps was printed to a file in the job dir
err_path = paste(job_dir,"chr_filter.err",sep="")
log_path = paste(job_dir,"chr_filter.log",sep="")
chr_filter = "--chr 0-22"
if(!autosomal_chrs){
  chr_filter = "--no-chr 0-22"
}
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data",sep=''),
                 chr_filter,
                 "--missing --freq --make-bed --out",paste(job_dir,"merged_mega_data_autosomal",sep=''))
curr_sh_file = "chr_filter.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(120)
print("After removing non autosomal variants, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.bim",sep="")))))
ids = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]!="Cooper")))
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Cooper",na.rm = T)))
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Elite",na.rm = T)))

# quick comparison of the data using TOP strand vs + strand
# plus_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_fwd_strand/raw.bim")
# top_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl/raw.bim")
# dim(plus_bim)
# dim(top_bim)
# comparisons = list()
# for(j in 1:ncol(plus_bim)){
#   comparisons[[j]] = table(plus_bim[,j]==top_bim[,j])  
# }

####################################################################################################
####################################################################################################
####################################################################################################
# Here we analyze the results of some of the jobs above to determine the final set of subjects
# for the analysis. We basically exclude sex check failures and low call rate samples.

# compare to known sex from the metadata
metadata_sex = sample_metadata_raw$Sex_input_data
names(metadata_sex) = rownames(sample_metadata_raw)

imputed_sex1 = read_plink_table(paste(input_bfile1,".sexcheck",sep=''))[,4]
imputed_sex2 = read_plink_table(paste(input_bfile2,".sexcheck",sep=''))[,4]
inds1 = intersect(names(imputed_sex1),names(metadata_sex))
inds2 = intersect(names(imputed_sex2),names(metadata_sex))

sex_errs1 = inds1[imputed_sex1[inds1] == "0" | ((imputed_sex1[inds1] !="2" & metadata_sex[inds1]=="F") | (imputed_sex1[inds1] !="1"& metadata_sex[inds1]=="M"))]
sex_errs2 = inds2[imputed_sex2[inds2] == "0" | ((imputed_sex2[inds2] !="2" & metadata_sex[inds2]=="F") | (imputed_sex2[inds2] !="1"& metadata_sex[inds2]=="M"))]
sex_errs = union(sex_errs1,sex_errs2)
ids = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("sex errors on the raw data, number of errors:",length(sex_errs)))
sex_errs = intersect(sex_errs,ids)
print(paste("number of errors in samples that survived prev filters:",length(sex_errs)))

# Check the ids and compare to the exome data
# Create a report with the sex errors
m = sample_metadata_raw[sex_errs,]$Cohort
m = cbind(sample_metadata_raw[sex_errs,"Sample_ID"],m)
rownames(m) = sex_errs
colnames(m) = c("Sample_ID","Cohort")
write.table(m,file=paste(job_dir,"sex_impute_analysis_report.txt",sep=''),sep="\t",quote=F)

# Missigness report after quality and maf filtering
# look at the results, compare to Illumina's
missinigness_report = read_plink_table(paste(job_dir,"merged_mega_data_autosomal.imiss",sep=''))
call_rates_after_filters = 1-as.numeric(missinigness_report[,6])
names(call_rates_after_filters) = rownames(missinigness_report)
low_cr_samples = missinigness_report[call_rates_after_filters<final_subj_min_call_rate,2]
length(intersect(low_cr_samples,sex_errs))

curr_fam = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)
final_subject_qc_excluded_samples = union(sex_errs,low_cr_samples)
curr_fam = curr_fam[is.element(curr_fam[,2],set=final_subject_qc_excluded_samples),]
remove_subjects_using_plink(paste(job_dir,"merged_mega_data_autosomal",sep=""),
                            curr_fam,
                            job_dir,"final_subj_qc","merged_mega_data_autosomal",
                            batch_script_func=get_sh_default_prefix)
wait_for_job(60)
# Report low call rate samples
m = sample_metadata_raw[low_cr_samples,]$Cohort
# look at call rates of the samples with error in sex imputation
m = cbind(m,call_rates_after_filters[low_cr_samples])
m = cbind(sample_metadata_raw[low_cr_samples,"Sample_ID"],m)
colnames(m) = c("Sample_ID","Cohort","Call_rate_qc_maf_filters")
write.table(m,file=paste(job_dir,"missigness_analysis_report.txt",sep=''),sep="\t",quote=F)

print("After removing non autosomal variants, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.bim",sep="")))))
ids = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Cooper",na.rm = T)))
print(paste("number of ELITE samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="ELITE",na.rm = T)))
print(paste("number of GP samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="genepool",na.rm = T)))

# Organize the new phenotypes

# Add imputed sex
imputed_sex1 = read_plink_table(paste(input_bfile1,".sexcheck",sep=''))[,4]
imputed_sex2 = read_plink_table(paste(input_bfile2,".sexcheck",sep=''))[,4]
fam_samples = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)
subjects_for_analysis = fam_samples[,2]
is_mega_consortium = !grepl(fam_samples[,1],pattern="^file")
imputed_sex = c(
  imputed_sex1[is.element(names(imputed_sex1),set=subjects_for_analysis[!is_mega_consortium])],
  imputed_sex2[is.element(names(imputed_sex2),set=subjects_for_analysis[is_mega_consortium])]
)
length(intersect(names(imputed_sex),subjects_for_analysis)) == length(subjects_for_analysis) # must be TRUE
imputed_sex = imputed_sex[subjects_for_analysis]
all(names(imputed_sex) == fam_samples[,2])

final_sample_set_pheno = cbind(sample_metadata_raw[ids,],missinigness_report[ids,"F_MISS"],imputed_sex[ids])
save(final_sample_set_pheno,file=paste(job_dir,"final_sample_set_pheno.RData",sep=""))

# Samples: put all filtered samples in one file 
load("final_sample_set_pheno.RData")
sample_metadata_raw_ids = paste(sample_metadata_raw[,1],sample_metadata_raw[,2],sep="_")
S1 = intersect(rownames(sample_data1),sample_metadata_raw_ids)
S2 = intersect(rownames(sample_data2),sample_metadata_raw_ids)
removed1 = setdiff(S1,rownames(final_sample_set_pheno))
removed2 = setdiff(S2,rownames(final_sample_set_pheno))
gp_samples = sample_metadata_raw_ids[sample_metadata_raw$Cohort=="genepool"]
removed1 = setdiff(removed1,gp_samples)
removed2 = setdiff(removed2,gp_samples)
removed = union(removed1,removed2)
xx = sample_metadata_raw
removed_info = xx[is.element(sample_metadata_raw_ids,removed),1:10]
write.table(removed_info,file="removed_samples_info.txt",sep="\t",row.names = F,col.names = T,quote=F)

####################################################################################################
####################################################################################################
####################################################################################################
# GWAS vs. sex
imputed_sex1 = read_plink_table(paste(input_bfile1,".sexcheck",sep=''))[,4]
imputed_sex2 = read_plink_table(paste(input_bfile2,".sexcheck",sep=''))[,4]
fam_samples = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)
subjects_for_analysis = fam_samples[,2]
is_mega_consortium = !grepl(fam_samples[,1],pattern="^file")
imputed_sex = c(
  imputed_sex1[is.element(names(imputed_sex1),set=subjects_for_analysis[!is_mega_consortium])],
  imputed_sex2[is.element(names(imputed_sex2),set=subjects_for_analysis[is_mega_consortium])]
)
length(intersect(names(imputed_sex),subjects_for_analysis)) == length(subjects_for_analysis) # must be TRUE
imputed_sex = imputed_sex[subjects_for_analysis]
all(names(imputed_sex) == fam_samples[,2])

pheno_file = paste(job_dir,"fam_with_sex_info.txt",sep='')
sex_fam_info = cbind(fam_samples[,1:2],imputed_sex)
colnames(sex_fam_info) = c("FID","IID","sex")
write.table(sex_fam_info,file=pheno_file,row.names = F,sep="\t",col.names = T,quote = F)
err_path = paste(job_dir,"sex_gwas_qc.err",sep="")
log_path = paste(job_dir,"sex_gwas_qc.log",sep="")
curr_cmd = paste("plink2",
                 "--bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
                 "--logistic hide-covar firth-fallback",
                 paste("--pheno",pheno_file),
                 paste("--pheno-name sex"),
                 "--adjust",
                 "--out",paste(job_dir,"sex_gwas_qc",sep=''))
curr_sh_file = "sex_gwas_qc.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(120)

# Look at the results 
sex_gwas_qc_res = read.table(paste(job_dir,"sex_gwas_qc.sex.glm.logistic.hybrid.adjusted",sep=""),stringsAsFactors = F)
table(sex_gwas_qc_res[,3]> 1e-5)
sex_snps = sex_gwas_qc_res[sex_gwas_qc_res[,3]> 1e-5,2]
extract_snps_using_plink(paste(job_dir,"merged_mega_data_autosomal",sep=""),sex_snps,job_dir,
                         "_snps_to_keep_after_sex_gwas_qc","merged_mega_data_autosomal",
                         batch_script_func=get_sh_default_prefix)
wait_for_job(120)

print("After removing non autosomal variants, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal.bim",sep="")))))
ids = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Cooper",na.rm = T)))
print(paste("number of ELITE samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="ELITE",na.rm = T)))
print(paste("number of GP samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="genepool",na.rm = T)))

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # GWAS vs. call rate
# pheno_file = paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep='')
# err_path = paste(job_dir,"call_rate_gwas_qc.err",sep="")
# log_path = paste(job_dir,"call_rate_gwas_qc.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
#                  "--linear hide-covar",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name call_rates_after_filters[subjects_for_analysis]"),
#                  "--adjust",
#                  "--out",paste(job_dir,"call_rate_gwas_qc",sep=''))
# curr_sh_file = "call_rate_gwas_qc.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# 
# # Look at the results 
# call_rate_gwas_qc_res = read.table(paste(job_dir,"call_rate_gwas_qc.call_rates_after_filters[subjects_for_analysis].glm.linear.adjusted",sep=""),stringsAsFactors = F)
# rownames(call_rate_gwas_qc_res) = call_rate_gwas_qc_res[,2]
# call_rate_gwas_qc_res["rs1747677",]
# call_rate_gwas_qc_res["rs73017651",]
# call_rate_gwas_qc_res["rs867306",]
# call_rate_gwas_qc_res["rs137925455",]
# # rs73017651 rs867306 rs137925455
# callrate_snps = call_rate_gwas_qc_res[call_rate_gwas_qc_res[,3]> 1e-4,2]
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Create a freq file for each cohort
# TODO:  move to a later stage
#
# cohorts = sample_metadata_raw$Cohort
# names(cohorts) = rownames(sample_metadata_raw)
# fam = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=''),stringsAsFactors = F)
# rownames(fam) = fam[,2]
# cohorts = cohorts[rownames(fam)]
# 
# for(cc in unique(cohorts)){
#   inds = cohorts == cc
#   m = fam[inds,1:2]
#   write.table(m,sep=" ",file=paste(job_dir,cc,"_subjects.txt",sep=""),row.names = F,quote = F)
#   err_path = paste(job_dir,cc,"_cohort_freq.err",sep="")
#   log_path = paste(job_dir,cc,"_cohort_freq.log",sep="")
#   curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
#                    "--keep",paste(job_dir,cc,"_subjects.txt",sep=""),
#                    "--freq --out",paste(job_dir,cc,"_cohort_freq",sep=""))
#   curr_sh_file = paste(cc,"_cohort_freq.sh",sep="")
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_default_prefix(err_path,log_path),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# }
# wait_for_job()
# 
# # analyze the results, get low MAF SNPs
# mafs = list()
# get_mafs<-function(path){
#   x = read.table(path,header = T,stringsAsFactors = F,row.names = 2)
#   y = x[,"MAF"]
#   names(y) = rownames(x)
#   return(y)
# }
# mafs[["all"]] = get_mafs(paste(job_dir,"merged_mega_data_autosomal.frq",sep=''))
# for(cc in unique(cohorts)){
#   mafs[[as.character(cc)]] = get_mafs(paste(job_dir,cc,"_cohort_freq.frq",sep=""))
#   print(all(names(mafs[[as.character(cc)]])==names(mafs[["all"]])))
# }
# mafs = sapply(mafs,function(x)x) # into a matrix
# table(apply(mafs <0.001,1,any))
# # > cor(mafs,method="spearman")
# # all         2         1
# # all 1.0000000 0.9955802 0.9981975
# # 2   0.9955802 1.0000000 0.9884692
# # 1   0.9981975 0.9884692 1.0000000
# mafs_0.5 = mafs < 0.05
# table(rowSums(mafs_0.5))
# table(mafs_0.5[,1],mafs_0.5[,2] | mafs_0.5[,3])
# 
# high_maf_snps_for_analysis = rownames(mafs)[mafs[,2]>= 0.05 & mafs[,3]>= 0.05]
# extract_snps_using_plink(paste(job_dir,"merged_mega_data_autosomal",sep=""),high_maf_snps_for_analysis,job_dir,
#                          "high_maf_snps_for_analysis","merged_mega_data_autosomal_after_maf",
#                          batch_script_func=get_sh_default_prefix)

####################################################################################################
####################################################################################################
####################################################################################################
# Prune and run PCA and relatedness
analysis_name = "our_data_ld_prune"
err_path = paste(job_dir,analysis_name,"_ld_report.err",sep="")
log_path = paste(job_dir,analysis_name,"_ld_report.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
                 "--maf 0.05 --indep-pairwise 250 10",0.1,
                 "--out",paste(job_dir,analysis_name,"_plink.prune",sep=""))
curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(60)

print(paste("Prune before PCA, num of variants is:",
            length(readLines(paste(job_dir,analysis_name,"_plink.prune.prune.in",sep="")))))

err_path = paste(job_dir,"final_data_pca.err",sep="")
log_path = paste(job_dir,"final_data_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
                 "--extract", paste(job_dir,analysis_name,"_plink.prune.prune.in",sep=""),
                 "--pca 40 --out",paste(job_dir,"merged_mega_data_autosomal",sep=''))
curr_sh_file = "final_data_pca.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# relatedness
err_path = paste(job_dir,"final_data_rl.err",sep="")
log_path = paste(job_dir,"final_data_rl.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
                 "--extract", paste(job_dir,analysis_name,"_plink.prune.prune.in",sep=""),
                 "--genome --min 0.2 --out",paste(job_dir,"merged_mega_data_autosomal",sep=''))
curr_sh_file = "final_data_rl.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# frq on all
err_path = paste(job_dir,"final_data_fr.err",sep="")
log_path = paste(job_dir,"final_data_fr.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal",sep=''),
                 "--freq --out",paste(job_dir,"merged_mega_data_autosomal",sep=''))
curr_sh_file = "final_data_fr.sh"
print_sh_file(paste(job_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
wait_for_job(120)

# Code below wraps all relevant metadata and print it into a file
imputed_sex1 = read_plink_table(paste(input_bfile1,".sexcheck",sep=''))[,4]
imputed_sex2 = read_plink_table(paste(input_bfile2,".sexcheck",sep=''))[,4]
fam_samples = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)
pca_res = read_pca_res(paste(job_dir,"merged_mega_data_autosomal.eigenvec",sep=""))
all(is.element(fam_samples[,2],set=rownames(sample_metadata_raw))) # this should be TRUE
all(rownames(pca_res) == fam_samples[,2]) # this should be true as well
subjects_for_analysis = fam_samples[,2]
is_mega_consortium = !grepl(fam_samples[,1],pattern="^file")
imputed_sex = c(
  imputed_sex1[is.element(names(imputed_sex1),set=subjects_for_analysis[!is_mega_consortium])],
  imputed_sex2[is.element(names(imputed_sex2),set=subjects_for_analysis[is_mega_consortium])]
)
length(intersect(names(imputed_sex),subjects_for_analysis)) == length(subjects_for_analysis) # must be TRUE
imputed_sex = imputed_sex[subjects_for_analysis]

missinigness_report = read_plink_table(paste(job_dir,"merged_mega_data_autosomal.imiss",sep=''))
call_rates_after_filters = 1-as.numeric(missinigness_report[,6])
names(call_rates_after_filters) = rownames(missinigness_report)

# put all covariates in one table
covariate_matrix = cbind(fam_samples[,1:2],sample_metadata_raw[subjects_for_analysis,],
                         call_rates_after_filters[subjects_for_analysis],imputed_sex,is_mega_consortium,
                         pca_res)

# Correct some columns to make them easier to work with plink
covariate_matrix[covariate_matrix==""] = NA
covariate_matrix[covariate_matrix=="_"] = NA
covariate_matrix$Shipment.date[is.na(covariate_matrix$Shipment.date)] = "uknown"
covariate_matrix$Cohort[covariate_matrix$Cohort=="ELITE"] = "2"
covariate_matrix$Cohort[covariate_matrix$Cohort=="Cooper"] = "1"
covariate_matrix[,"Shipment.date"] = cov_phe_col_to_plink_numeric_format(covariate_matrix[,"Shipment.date"])
colnames(covariate_matrix)[colnames(covariate_matrix)=="Age..at.test."] = "age"
colnames(covariate_matrix)[colnames(covariate_matrix)=="imputed_sex"] = "sex"
colnames(covariate_matrix)[colnames(covariate_matrix)=="Shipment.date"] = "batch"
colnames(covariate_matrix)[1:2] = c("FID","IID")
covariate_matrix$batch = paste("batch",covariate_matrix$batch,sep="")
for(j in 1:ncol(covariate_matrix)){
  if(class(covariate_matrix[[j]]) == "character"){
    covariate_matrix[[j]] = gsub(pattern = " ",replacement = "_",x = covariate_matrix[[j]])
  }
}

write.table(covariate_matrix,file=
              paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),
            sep="\t",quote=F,row.names = F)
write.table(covariate_matrix,file=
              paste(job_dir,"integrated_sample_metadata_and_covariates.phe",sep=''),
            sep=" ",quote=F,row.names = F)

# Some stats
table(covariate_matrix[,"Cohort"])
table(covariate_matrix[,"Cohort"],covariate_matrix[,"batch"])

# !! At this point we have our directly genotyped data. Subsequent analyses can follow but this
# is the data we should use for the initial QCs, especially when integrating with external datasets.

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
load(paste(job_dir,"manual_clustering.RData",sep=""))
selected_subjects = names(which(manual_clustering))
curr_fam = read.table(paste(job_dir,"merged_mega_data_autosomal.fam",sep=""),stringsAsFactors = F)
curr_fam = curr_fam[!is.element(curr_fam[,2],set=selected_subjects),]
remove_subjects_using_plink(paste(job_dir,"merged_mega_data_autosomal",sep=""),
                            curr_fam,
                            job_dir,"_eu_selection","merged_mega_data_autosomal_eu_selected",
                            batch_script_func=get_sh_default_prefix)
wait_for_job(120)
print("After selecting the EUs, data sizes are:")
print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_eu_selected.fam",sep="")))))
print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_eu_selected.bim",sep="")))))
ids = read.table(paste(job_dir,"merged_mega_data_autosomal_eu_selected.fam",sep=""),stringsAsFactors = F)[,2]
print(paste("number of ELITE samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="ELITE",na.rm = T)))
print(paste("number of cooper samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="Cooper",na.rm = T)))
print(paste("number of genepool samples in this file:",sum(sample_metadata_raw[ids,"Cohort"]=="genepool",na.rm = T)))

curr_gp_data = sample_metadata_raw[ids,]
curr_gp_data = curr_gp_data[curr_gp_data[,"Cohort"]=="genepool",]
table(is.na(curr_gp_data$Age..at.test.)) 

# Load covariates, create a new dir for GWAS and reruns of PCA
curr_dir = paste(job_dir,"eu_gwas/",sep="")
system(paste("mkdir",curr_dir))
analysis_name = "our_data_ld_prune"
err_path = paste(curr_dir,analysis_name,"_ld_report.err",sep="")
log_path = paste(curr_dir,analysis_name,"_ld_report.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                 "--maf 0.05 --indep-pairwise 250 10",0.1,
                 "--out",paste(curr_dir,analysis_name,"_plink.prune",sep=""))
curr_sh_file = paste(analysis_name,"_ld_report.sh",sep="")
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
wait_for_job(60)

print(paste("Prune before PCA, num of variants is:",
            length(readLines(paste(curr_dir,analysis_name,"_plink.prune.prune.in",sep="")))))

err_path = paste(curr_dir,"final_data_pca.err",sep="")
log_path = paste(curr_dir,"final_data_pca.log",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                 "--extract", paste(curr_dir,analysis_name,"_plink.prune.prune.in",sep=""),
                 "--pca 40 --out",paste(curr_dir,"merged_mega_data_autosomal",sep=''))
curr_sh_file = "final_data_pca.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

# Update the covariates
covariate_matrix = read.table(paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),
                              sep="\t",header=T,stringsAsFactors = F)
rownames(covariate_matrix) = covariate_matrix[,2]
covariate_matrix = covariate_matrix[ids,]
cohort = covariate_matrix$Cohort
elite_vs_gp = cohort
elite_vs_gp[elite_vs_gp=="1"] = NA
elite_vs_gp[elite_vs_gp=="genepool"] = 1
elite_vs_gp[elite_vs_gp=="2"] = 2
table(elite_vs_gp)
cooper_vs_gp = cohort
cooper_vs_gp[cooper_vs_gp=="2"] = NA
cooper_vs_gp[cooper_vs_gp=="1"] = 2
cooper_vs_gp[cooper_vs_gp=="genepool"] = 1
table(cooper_vs_gp)
covariate_matrix = cbind(covariate_matrix,elite_vs_gp,cooper_vs_gp)

new_pca = read_pca_res(paste(curr_dir,"merged_mega_data_autosomal.eigenvec",sep=""))
# sanity check
nrow(new_pca) == length(intersect(rownames(new_pca),rownames(covariate_matrix)))
new_pca = new_pca[rownames(covariate_matrix),]
covariate_matrix[,paste("PC",1:40,sep="")] = new_pca
covariate_matrix[covariate_matrix=="_"] = NA
write.table(covariate_matrix,file=paste(curr_dir,"all_covs_and_pheno.phe",sep=""),row.names=F,
                                        col.names = T,sep=" ",quote=F)

# Run the GWASs
covs_file = paste(curr_dir,"all_covs_and_pheno.phe",sep="")
curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                 "--logistic hide-covar",
                 "--pheno",covs_file,
                 "--pheno-name elite_vs_gp",
                 "--covar",covs_file,
                 "--maf 0.01",
                 "--covar-name sex,age,PC1,PC2,PC3,PC4,PC5",
                 "--allow-no-sex --adjust",
                 "--threads",8,
                 "--out",paste(curr_dir,"elite_vs_gp_gwas_res",sep="")
)
run_plink_command(curr_cmd,curr_dir,paste("elite_vs_gp_gwas_res",sep=""),
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)

curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep=''),
                 "--logistic hide-covar",
                 "--pheno",covs_file,
                 "--pheno-name cooper_vs_gp",
                 "--covar",covs_file,
                 "--maf 0.01",
                 "--covar-name sex,age,PC1,PC2,PC3,PC4,PC5",
                 "--allow-no-sex --adjust",
                 "--threads",8,
                 "--out",paste(curr_dir,"cooper_vs_gp_gwas_res",sep="")
)
run_plink_command(curr_cmd,curr_dir,paste("cooper_vs_gp_gwas_res",sep=""),
                  get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=8,mem_size=32000)

# Print in fuma-compatible format
out_files = list.files(curr_dir)
out_files = out_files[grepl("logistic$",out_files)]
for(ff in out_files){
  res = read.table(paste(curr_dir,ff,sep=""),stringsAsFactors = F,header=T)
  print("####")
  print(table(as.numeric(res[,"P"]<1e-7)))
  print(table(as.numeric(res[,"P"]<5e-8)))
  print("####")
  fuma_res = res[,c("CHR","BP","P")]
  colnames(fuma_res) = c("chromosome","position","P-value")
  write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
              sep=" ")
}

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Take the mega direct geno as is, EU samples and run GWAS for each PC
# curr_dir = paste(job_dir,"eu_gwas/",sep="")
# covs_file = paste(curr_dir,"all_covs_and_pheno.phe",sep="")
# curr_bfile = paste(job_dir,"merged_mega_data_autosomal_eu_selected",sep='')
# for (j in 1:10){
#   currname = paste("PC",j,sep="")
#   curr_cmd = paste("plink --bfile",curr_bfile,
#                    "--linear hide-covar",
#                    "--pheno",covs_file,
#                    "--pheno-name", paste("PC",j,sep=""),
#                    "--covar",covs_file,
#                    "--maf 0.01",
#                    "--covar-name sex",
#                    "--allow-no-sex --adjust",
#                    "--threads",4,
#                    "--out",paste(curr_dir,currname,sep="")
#   )
#   run_plink_command(curr_cmd,curr_dir,currname,
#                     get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=16000)
# }
# 
# # Create FUMA files
# out_files = list.files(curr_dir)
# out_files = out_files[grepl("linear$",out_files)]
# out_files = out_files[!grepl("fuma",out_files)]
# for(ff in out_files){
#   res = read.table(paste(curr_dir,ff,sep=""),stringsAsFactors = F,header=T)
#   fuma_res = res[,c("CHR","BP","P")]
#   colnames(fuma_res) = c("chromosome","position","P-value")
#   write.table(fuma_res,file=paste(curr_dir,"fuma_",ff,sep=""),quote=F,col.names = T,row.names = F,
#               sep=" ")
#   print(ff)
# }
# 
# # Check the association with the cohort
# d = read.table(covs_file,header=T)
# pc_ps = c()
# for(j in 1:40){
#   curr_inds = d[,"Cohort"] == "2" | d[,"Cohort"]=="1"
#   p1 = compute_pc_vs_binary_variable_association_p(
#     pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"],
#     test=wilcox.test
#   )
#   curr_inds = d[,"Cohort"] == "2" | d[,"Cohort"]=="genepool"
#   p2 = compute_pc_vs_binary_variable_association_p(
#     pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"],
#     test=wilcox.test
#   )
#   curr_inds = d[,"Cohort"] == "1" | d[,"Cohort"]=="genepool"
#   p3 = compute_pc_vs_binary_variable_association_p(
#     pc = d[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"],
#     test=wilcox.test
#   )
#   pc_ps = rbind(pc_ps,c(p1,p2,p3))
# }
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################

####################################################################################################
####################################################################################################
####################################################################################################
# Transform the dataset into HRC-based or 1000G-based data
# Run the check_bim analysis

curr_dir = paste(job_dir,"1000g/",sep="")
system(paste("mkdir",curr_dir))
setwd(curr_dir)
curr_bfile = "merged_mega_data_autosomal"
check_if_bim_is_sorted(paste(paste(job_dir,curr_bfile,".bim",sep='')))
err_path = paste(curr_dir,"run_check_bim.err",sep="")
log_path = paste(curr_dir,"run_check_bim.log",sep="")
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",curr_dir))
# For 1000G-based analysis
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(job_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(job_dir,curr_bfile,".frq",sep=''),
                 "-1000g -p EUR -t 0.3 -r ",
                 "/home/users/davidama/apps/check_bim/1000GP_Phase3_combined.legend")
curr_sh_file = "run_check_bim.sh"
# print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
#               get_sh_prefix_bigmem(err_path,log_path,mem_size = 256000,time="6:00:00"),curr_cmd)
# Try without bigmem
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
wait_for_job()
wait_for_job()
system(paste("less ",curr_dir,"Run-plink.sh | grep TEMP > ",curr_dir,"Run-plink_1000g.sh",sep=""))
run_sh_lines = readLines(paste(curr_dir,"Run-plink_1000g.sh",sep=""))
run_sh_lines[1] = gsub(paste("plink --bfile",curr_bfile),paste("plink --bfile ",job_dir,curr_bfile,sep=""),run_sh_lines[1])
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "")
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

# For HRC-based analysis
curr_dir = paste(job_dir,"hrc/",sep="")
system(paste("mkdir",curr_dir))
setwd(curr_dir)
system(paste("cp /home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl",curr_dir))
err_path = paste(curr_dir,"run_check_bim2.err",sep="")
log_path = paste(curr_dir,"run_check_bim2.log",sep="")
curr_cmd = paste("perl", paste(curr_dir, "HRC-1000G-check-bim-NoReadKey.pl",sep=""),
                 "-b", paste(job_dir,curr_bfile,".bim",sep=''),
                 "-f", paste(job_dir,curr_bfile,".frq",sep=''),
                 "-hrc -p EU -t 0.3 -r",
                 "/home/users/davidama/apps/check_bim/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
curr_sh_file = "run_check_bim2.sh"
# print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
#               get_sh_prefix_bigmem(err_path,log_path,mem_size = 256000,time="6:00:00"),curr_cmd)
# Try without bigmem
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,
              mem_size = 64000,time="6:00:00",Ncpu = 4),curr_cmd)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))
wait_for_job()
system(paste("less ",curr_dir,"Run-plink.sh | grep TEMP > ",curr_dir,"Run-plink_hrc.sh",sep=""))
run_sh_lines = readLines(paste(curr_dir,"Run-plink_hrc.sh",sep=""))
run_sh_lines[1] = gsub(paste("plink --bfile",curr_bfile),paste("plink --bfile ",job_dir,curr_bfile,sep=""),run_sh_lines[1])
run_sh_lines = sapply(run_sh_lines,gsub,pattern = "-updated",replacement = "")
err_path = paste(curr_dir,"run_check_bim_update.err",sep="")
log_path = paste(curr_dir,"run_check_bim_update.log",sep="")
curr_sh_file = "run_check_bim_update.sh"
print_sh_file(paste(curr_dir,curr_sh_file,sep=''),
              get_sh_default_prefix(err_path,log_path),run_sh_lines)
system(paste("sbatch",paste(curr_dir,curr_sh_file,sep='')))

# To download and install the tools on the cluster
# 1. Check bim:
# wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# unzip HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# mkdir check_bim
# mv HRC* check_bim/
# mv LICENSE.txt check_bim/
# cd check_bim
# wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
# gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
#
# 2. Strand analysis:
# mkdir ~/apps/wrayner_strand/
# cd ~/apps/wrayner_strand
# wget http://www.well.ox.ac.uk/~wrayner/strand/update_build.sh

####################################################################################################
####################################################################################################
####################################################################################################

# !!!!!!!!!!!!!!!!!!!!
# Code beyond this point should be moved to another script
# We stop here because at this point we have MEGA EU data after standard QC steps

# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Run GWAS for each PC
# pheno_file = paste(job_dir,"integrated_sample_metadata_and_covariates.phe",sep='')
# curr_bfile = "merged_mega_data_autosomal_after_maf"
# for (j in 1:5){
#   err_path = paste(job_dir,"all_mega_data_gwas_PC",j,".err",sep="")
#   log_path = paste(job_dir,"all_mega_data_gwas_PC",j,".log",sep="")
#   curr_cmd = paste("plink2",
#                    "--bfile",paste(job_dir,curr_bfile,sep=''),
#                    "--linear hide-covar",
#                    paste("--pheno-name",paste("PC",j,sep="")),
#                    paste("--pheno",pheno_file),
#                    paste("--covar",pheno_file),
#                    "--covar-name sex,age",
#                    "--adjust",
#                    "--out",paste(job_dir,"all_mega_data_gwas_PC",j,"",sep=''))
#   curr_sh_file = paste("all_mega_data_gwas_PC",j,".sh",sep="")
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# }
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Use PCA and relatedness to remove samples and then rerun PCA
# covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates.txt",sep=''),stringsAsFactors = F)
# rownames(covariate_matrix) = covariate_matrix$IID
# altsamp_id = sample_metadata_raw$alt_sample_id
# names(altsamp_id) = rownames(sample_metadata_raw)
# is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = rownames(sample_metadata_raw)
# 
# set.seed(123)
# d = covariate_matrix
# pc_x = as.matrix(d[,paste("PC",1:3,sep="")])
# rownames(pc_x) = rownames(d)
# # kmeans_res = kmeans(pc_x,4)$cluster
# # Using hierarchical
# dd = dist(pc_x,method="manhattan")
# h = hclust(dd,method = "single")
# kmeans_res = run_hclust(pc_x,12,dd,h)
# 
# table(kmeans_res)
# table(kmeans_res,d$Cohort)
# table(kmeans_res,is_jap[names(kmeans_res)])
# 
# to_rem = rep(F,nrow(d))
# for(j in 1:20){
#   x = d[,paste("PC",j,sep="")]
#   x = (x-mean(x))/sd(x)
#   print(sum(abs(x)>6))
#   to_rem[abs(x)>6] = T
# }
# table(to_rem)
# 
# # Select subjects from the largest cluster
# clustable = table(kmeans_res)
# selected_cluster = names(which(clustable == max(clustable)))
# selected_subjects = names(which(kmeans_res == selected_cluster))
# 
# print(paste("PCA and clustering analysis, largest cluster size is:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# selected_subjects = setdiff(selected_subjects,rownames(d)[to_rem])
# 
# print(paste("PC outliers removed, number of remaining subjects:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# library("igraph",lib.loc = "~/R/packages")
# rl_data = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf.genome",sep=""),
#                      header=T,stringsAsFactors = F)
# rl_edges = as.matrix(rl_data[,c("IID1","IID2")])
# rl_g = igraph::graph_from_edgelist(rl_edges,directed = F)
# rl_clusters = clusters(rl_g)[[1]]
# rl_subjects_to_remove = c()
# for(cl in unique(rl_clusters)){
#   curr_subjects = names(rl_clusters)[rl_clusters==cl]
#   rl_subjects_to_remove = c(rl_subjects_to_remove,curr_subjects[-1])
# }
# print(paste("Relatedness analysis, number of subjects to remove:",length(rl_subjects_to_remove)))
# print(paste("Intersection with largest cluster",
#             length(intersect(rl_subjects_to_remove,selected_subjects))))
# selected_subjects = setdiff(selected_subjects,rl_subjects_to_remove)
# print(paste("PCA and clustering analysis, largest cluster size is:",length(selected_subjects)))
# print(paste("number of cooper subjects",sum(covariate_matrix[selected_subjects,"Cohort"]==1)))
# 
# curr_fam = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf.fam",sep=""),stringsAsFactors = F)
# curr_fam = curr_fam[!is.element(curr_fam[,2],set=selected_subjects),]
# remove_subjects_using_plink(paste(job_dir,"merged_mega_data_autosomal_after_maf",sep=""),
#                             curr_fam,
#                             job_dir,"_pca_and_rl_subj_qc","merged_mega_data_autosomal_after_maf_after_pca",
#                             batch_script_func=get_sh_default_prefix)
# wait_for_job()
# print("After PCA and Rl analysis, data sizes are:")
# print(paste("number of samples:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.fam",sep="")))))
# print(paste("number of snps:",length(readLines(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.bim",sep="")))))
# 
# # Prune and run PCA and relatedness
# analysis_name = "our_data_ld_prune2"
# err_path = paste(job_dir,analysis_name,"_ld_report2.err",sep="")
# log_path = paste(job_dir,analysis_name,"_ld_report2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--indep-pairwise 250 10",0.1,
#                  "--out",paste(job_dir,analysis_name,"_after_pca",sep=""))
# curr_sh_file = paste(analysis_name,"_ld_report2.sh",sep="")
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,Ncpu=2,mem_size=16000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# print(paste("Prune before PCA, num of variants is:",
#             length(readLines(paste(job_dir,analysis_name,"_after_pca.prune.in",sep="")))))
# err_path = paste(job_dir,"final_data_pca2.err",sep="")
# log_path = paste(job_dir,"final_data_pca2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--extract", paste(job_dir,analysis_name,"_after_pca.prune.in",sep=""),
#                  "--pca --out",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''))
# curr_sh_file = "final_data_pca2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# # rerun frq
# err_path = paste(job_dir,"final_data_frq2.err",sep="")
# log_path = paste(job_dir,"final_data_frq2.log",sep="")
# curr_cmd = paste("plink --bfile",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''),
#                  "--freq --out",paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep=''))
# curr_sh_file = "final_data_frq2.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_default_prefix(err_path,log_path),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# 
# # Create the new covariate matrix file
# new_pca_res = read_pca_res(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.eigenvec",sep=""))
# covariate_matrix = covariate_matrix[rownames(new_pca_res),]
# covariate_matrix[,paste("PC",1:20,sep="")] = new_pca_res
# write.table(covariate_matrix,file=
#               paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.txt",sep=''),
#             sep="\t",quote=F,row.names = F)
# write.table(covariate_matrix,file=
#               paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.phe",sep=''),
#             sep=" ",quote=F,row.names = F)
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Quick GWAS runs between elite and genepool
# covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.txt",sep=''),
#                               stringsAsFactors = F)
# rownames(covariate_matrix) = covariate_matrix$IID
# 
# # read our fam file
# fam_info = read_plink_table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.fam",sep=""),has_header = F)
# iid_to_fid = fam_info[,1]
# 
# # PC vs cohort p-values
# pc_ps = c()
# for(j in 1:20){
#   pc_ps[j] = compute_pc_vs_binary_variable_association_p(
#     covariate_matrix[,paste("PC",j,sep="")],
#     covariate_matrix[,"Cohort"]
#   )
# }
# pc_ps = p.adjust(pc_ps)
# pc_ind = max(which(pc_ps<0.01))
# 
# # Logistic: Cooper vs. Elite: without batch but with the selected PCs
# pheno_file = paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.phe",sep='')
# curr_bfile = "merged_mega_data_autosomal_after_maf_after_pca"
# 
# err_path = paste(job_dir,"cooper_vs_elite.err",sep="")
# log_path = paste(job_dir,"cooper_vs_elite.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,curr_bfile,sep=''),
#                  "--logistic hide-covar firth-fallback",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name Cohort"),
#                  paste("--covar",pheno_file),
#                  paste("--covar-name sex,age,",paste(paste("PC",1:pc_ind,sep=""),collapse=","),sep=""),
#                  "--adjust",
#                  "--out",paste(job_dir,"cooper_vs_elite",sep=''))
# curr_sh_file = "cooper_vs_elite.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# 
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Run GWAS for each PC
# pheno_file = paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.phe",sep='')
# curr_bfile = "merged_mega_data_autosomal_after_maf_after_pca"
# for (j in 1:20){
#   err_path = paste(job_dir,"gwas_PC",j,".err",sep="")
#   log_path = paste(job_dir,"gwas_PC",j,".log",sep="")
#   curr_cmd = paste("plink2",
#                    "--bfile",paste(job_dir,curr_bfile,sep=''),
#                    "--linear hide-covar",
#                    paste("--pheno-name",paste("PC",j,sep="")),
#                    paste("--pheno",pheno_file),
#                    paste("--covar",pheno_file),
#                    "--covar-name sex,age",
#                    "--adjust",
#                    "--out",paste(job_dir,"gwas_PC",j,"",sep=''))
#   curr_sh_file = paste("gwas_PC",j,".sh",sep="")
#   print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#                 get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
#   system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# }
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Run cooper running times GWAS
# cooper_metadata = read.delim("/oak/stanford/groups/euan/projects/fitness_genetics/metadata/cooper_metadata.txt",
#                              stringsAsFactors = F,sep="\t")
# covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.txt",sep=''),
#                               stringsAsFactors = F)
# 
# # Create the pheno file for the analysis
# # 1. Parse the treadmill times
# xt = cooper_metadata$Treadmill.time
# new_xt = c()
# for(x in xt){
#   arr = strsplit(x,split=":")[[1]]
#   if(length(arr)>2){arr = arr[-3]}
#   new_xt = c(new_xt,
#              as.numeric(arr[1])+as.numeric(arr[2])/60)
# }
# names(new_xt) = as.character(cooper_metadata[,1])
# 
# new_xt = qnorm((rank(new_xt)-0.5)/length(new_xt))
# 
# # 2. Cut the cov matrix to have cooper samples only
# inds = is.element(covariate_matrix$Sample_ID,set=as.character(cooper_metadata[,1]))
# covariate_matrix = covariate_matrix[inds,]
# covariate_matrix[,"Treadmill.time"] = new_xt[as.character(covariate_matrix$Sample_ID)]
# table(is.na(covariate_matrix[,"Treadmill.time"]))
# 
# # 3. PC vs. time correlations
# pc_ps = c()
# for(j in 1:20){
#   x1 = covariate_matrix[,paste("PC",j,sep="")]
#   x2 = covariate_matrix[,"Treadmill.time"]
#   pc_ps[j] = cor.test(x1,x2,method="spearman")$p.value
# }
# 
# write.table(covariate_matrix,file=
#               paste(job_dir,"cooper_metadata_and_covariates.phe",sep=''),
#             sep=" ",quote=F,row.names = F)
# 
# # Run the GWAS
# pheno_file = paste(job_dir,"cooper_metadata_and_covariates.phe",sep='')
# err_path = paste(job_dir,"cooper_treadmill_times.err",sep="")
# log_path = paste(job_dir,"cooper_treadmill_times.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,curr_bfile,sep=''),
#                  "--linear hide-covar",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name Treadmill.time"),
#                  paste("--covar",pheno_file),
#                  paste("--covar-name sex,age,batch,",paste(paste("PC",1:3,sep=""),collapse=","),sep=""),
#                  "--adjust",
#                  "--out",paste(job_dir,"cooper_treadmill_times",sep=''))
# curr_sh_file = "cooper_treadmill_times.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # Run ELITE VO2 GWAS
# elite_metadata = read.delim("/oak/stanford/groups/euan/projects/fitness_genetics/metadata/elite_metadata.txt",
#                              stringsAsFactors = F,sep="\t")
# covariate_matrix = read.delim(paste(job_dir,"integrated_sample_metadata_and_covariates_after_pca1.txt",sep='')
#                               ,stringsAsFactors = F)
# 
# # change the vo2 feature names
# colnames(elite_metadata)[grepl("VO2",colnames(elite_metadata),ignore.case = T)] = c("vo2","vo2_l")
# new_xt = elite_metadata[,"vo2"]
# new_xt = qnorm((rank(new_xt)-0.5)/length(new_xt))
# names(new_xt) = elite_metadata[,"Sample_ID"]
# 
# # Cut the cov matrix to have cooper samples only
# inds = is.element(covariate_matrix$Sample_ID,set=as.character(elite_metadata[,1]))
# covariate_matrix = covariate_matrix[inds,]
# covariate_matrix[,"VO2max..ml.kg.min."] = new_xt[as.character(covariate_matrix$Sample_ID)]
# colnames(covariate_matrix)[grepl("VO2",colnames(covariate_matrix),ignore.case = T)] = c("vo2","vo2_l")
# covariate_matrix = covariate_matrix[!is.na(covariate_matrix[,"vo2"]),]
# quantile(covariate_matrix[,"vo2"])
# 
# # 3. PC vs. time correlations
# pc_ps = c()
# for(j in 1:20){
#   x1 = covariate_matrix[,paste("PC",j,sep="")]
#   x2 = covariate_matrix[,"vo2"]
#   pc_ps[j] = cor.test(x1,x2,method="spearman")$p.value
# }
# 
# write.table(covariate_matrix,file=
#               paste(job_dir,"elite_metadata_and_covariates.phe",sep=''),
#             sep=" ",quote=F,row.names = F)
# 
# # Run the GWAS
# pheno_file = paste(job_dir,"elite_metadata_and_covariates.phe",sep='')
# err_path = paste(job_dir,"elite_treadmill_times.err",sep="")
# log_path = paste(job_dir,"elite_treadmill_times.log",sep="")
# curr_cmd = paste("plink2",
#                  "--bfile",paste(job_dir,curr_bfile,sep=''),
#                  "--linear hide-covar",
#                  paste("--pheno",pheno_file),
#                  paste("--pheno-name vo2"),
#                  paste("--covar",pheno_file),
#                  paste("--covar-name sex,age,batch,",paste(paste("PC",1:3,sep=""),collapse=","),sep=""),
#                  "--adjust",
#                  "--out",paste(job_dir,"elite_treadmill_times",sep=''))
# curr_sh_file = "elite_treadmill_times.sh"
# print_sh_file(paste(job_dir,curr_sh_file,sep=''),
#               get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,"plink/2.0a1",2,10000),curr_cmd)
# system(paste("sbatch",paste(job_dir,curr_sh_file,sep='')))
# wait_for_job()
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # output all gwas results into a new dir with input files for fuma interpretation
# 
# res1 = read.table(paste(job_dir,"cooper_treadmill_times.Treadmill.time.glm.linear.adjusted",sep=""),
#                   stringsAsFactors = F)
# res2 = read.table(paste(job_dir,"elite_treadmill_times.vo2.glm.linear.adjusted",sep=""),
#                   stringsAsFactors = F)
# x1 = res1[,3];names(x1)=res1[,2]
# x2 = res2[,3];names(x2)=res2[,2]
# cor(x1,x2[names(x1)])
# fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
# new_ps = apply(cbind(x1,x2[names(x1)]),1,fishersMethod)
# 
# curr_bfile = paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca",sep='')
# create_fuma_files_for_fir(job_dir,
#                           paste(curr_bfile,".bim",sep=""),
#                           paste(curr_bfile,".frq",sep=""),p = 1,maf = 0.01,
#                           snps_to_exclude_from_results=NULL)
# 
# ####################################################################################################
# ####################################################################################################
# ####################################################################################################
# # QC steps: compare our datasets with other datasets
# 
# # for cohort specific mafs see _cohort_freq as regex
# 
# # 1. MAFs: ours vs. the ukbb
# our_mafs = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.frq",sep=""),
#                       header=T,stringsAsFactors = F)
# ukbb_mafs = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno.frq",
#                        header=T,stringsAsFactors = F)
# rownames(our_mafs) = our_mafs$SNP
# rownames(ukbb_mafs) = ukbb_mafs$SNP
# curr_snps = intersect(our_mafs$SNP,ukbb_mafs$SNP)
# x1 = our_mafs[curr_snps,"MAF"]
# x2 = ukbb_mafs[curr_snps,"MAF"]
# cor(x1,x2,method="spearman") # Spearman correlation: 0.997
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 441961    146 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 442034     73 
# 
# # 2. MAFs: ours vs. the ukbb after 1000g as ref
# our_mafs = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.frq",sep=""),
#                       header=T,stringsAsFactors = F)
# ukbb_mafs = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.frq",
#                        header=T,stringsAsFactors = F)
# rownames(our_mafs) = our_mafs$SNP
# rownames(ukbb_mafs) = ukbb_mafs$SNP
# curr_snps = intersect(our_mafs$SNP,ukbb_mafs$SNP)
# x1 = our_mafs[curr_snps,"MAF"]
# x2 = ukbb_mafs[curr_snps,"MAF"]
# cor(x1,x2,method="spearman") # Spearman correlation: 0.997
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 437340     85 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 437395     30
# y1 = our_mafs[curr_snps,c("A1","A2")]
# y1_flipped = t(apply(y1,1,flip_snp_info))
# y2 = ukbb_mafs[curr_snps,c("A1","A2")]
# # table(y1[,1]==y2[,1])
# # FALSE   TRUE 
# # 58274 379151 
# # table(y1[,2]==y2[,1])
# # FALSE   TRUE 
# # 432706   4719 
# # table(y1_flipped[,1]==y2[,1])
# # FALSE   TRUE 
# # 382415  55010
# strand_flipped_inds = y1_flipped[,1]==y2[,1]
# maf_flipped_inds = y1[,2]==y2[,1]
# # quantile(x2[maf_flipped_inds])
# # 0%      25%      50%      75%     100% 
# # 0.001459 0.170950 0.481500 0.494000 0.500000 
# # quantile(x2[strand_flipped_inds])
# # 0%       25%       50%       75%      100% 
# # 0.0009767 0.1382000 0.2552000 0.3751000 0.5000000 
# # quantile(x2)
# # 0%       25%       50%       75%      100% 
# # 0.0007508 0.0513200 0.1662000 0.3174000 0.5000000 
# 
# # Check the JHU stuff
# jhu_snps = our_mafs[grepl(our_mafs$SNP,pattern = "JHU"),]
# our_bim = read.table(paste(job_dir,"merged_mega_data_autosomal_after_maf_after_pca.bim",sep=""),
#                       header=F,stringsAsFactors = F)
# rownames(our_bim) = our_bim[,2]
# ukbb_bim = read.table("/oak/stanford/groups/euan/projects/fitness_genetics/ukbb/ukbb_imputed_20k_rand_controls_sex_age/merged_control_geno-1000g_updated.bim",
#                        header=F,stringsAsFactors = F)
# inds1 = is.element(ukbb_bim[,4],set = our_bim[jhu_snps$SNP,4])
# ukbb_bim = ukbb_bim[inds1,]
# jhu_bim = our_bim[jhu_snps$SNP,]
# rownames(ukbb_bim) = ukbb_bim[,2]
# ids1 = apply(jhu_bim[,c(1,4)],1,paste,collapse=";")
# ids2 = apply(ukbb_bim[,c(1,4)],1,paste,collapse=";")
# rownames(ukbb_bim) = ids2
# jhu_bim = jhu_bim[ids1!="0;0",]
# ids1 = ids1[ids1!="0;0"]
# rownames(jhu_bim) = ids1
# ids = intersect(ids1,ids2)
# jhu_bim = jhu_bim[ids,]
# ukbb_bim = ukbb_bim[ids,]
# x1 = our_mafs[jhu_bim[,2],"MAF"]
# x2 = ukbb_mafs[ukbb_bim[,2],"MAF"]
# # > table(abs(x1-x2)>0.1)
# # FALSE   TRUE 
# # 142586   1977 
# # > table(abs(x1-x2)>0.2)
# # FALSE   TRUE 
# # 143707    856 

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# Locally, should be commented out before running as a batch
setwd("/Users/David/Desktop/elite/november2018_analysis/mega_with_gp/")
d = read.delim("integrated_sample_metadata_and_covariates.txt",stringsAsFactors = F)
rownames(d) = d$IID
d2 = read.delim("../../metadata/june_2018_integrated_info/merged_metadata_file_stanford3k_elite_cooper.txt")
d2_ids = paste(d2$SentrixBarcode_A,d2$SentrixPosition_A,sep="_")
samp_id = d2$Sample_ID
altsamp_id = d2$alt_sample_id
names(samp_id) = d2_ids
names(altsamp_id) = d2_ids
is_jap = grepl(altsamp_id,pattern="JA"); names(is_jap) = d2_ids
table(d$Cohort)
d$Cohort[d$Cohort=="2"]="ELITE"
d$Cohort[d$Cohort=="1"] = "Cooper"

# Cluster by the first two PCs
set.seed(123)
pc_x = as.matrix(d[,paste("PC",1:2,sep="")])
rownames(pc_x) = rownames(d)

to_rem = rep(F,nrow(d))
for(j in 1:20){
  x = d[,paste("PC",j,sep="")]
  x = (x-mean(x))/sd(x)
  print(sum(abs(x)>6))
  to_rem[abs(x)>6] = T
}
table(to_rem)
# pc_x = pc_x[!to_rem,]

# Code for automatic clustering
# pc_x_kmeans = kmeans(pc_x,2)
# kmeans_res = pc_x_kmeans$cluster
# # Using hierarchical
# dd = dist(pc_x,method="euclidean")
# h = hclust(dd,method = "complete")
# kmeans_res = run_hclust(pc_x,5,dd,h)
# kmeans_res[kmeans_res!=1] = 0
# 
# kmeans_res[rownames(d)[to_rem]] = 10
# kmeans_res = kmeans_res[rownames(d)]
# table(kmeans_res)
# table(kmeans_res,d$Cohort)
# table(kmeans_res,is_jap[names(kmeans_res)])
# 
# res = two_d_plot_visualize_covariate(d$PC1[inds],
#   d$PC2[inds],kmeans_res[inds],kmeans_res[inds],
#   main = "PCA+Clustering",xlab="PC1",ylab="PC2",lwd=2,cex.axis=1.4,cex.lab=1.4,
#   xlim = c(-0.02,0.02),ylim = c(-0.05,0.1))
# legend(x="topright",names(res[[1]]),fill = res[[1]],cex=1.3,ncol = 4)


inds = 1:nrow(d)
res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC1",ylab="PC2",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="topleft",names(res[[1]]),
       fill = res[[1]],cex=1.1)

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC4[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC3",ylab="PC4",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="topright",names(res[[1]]),
       fill = res[[1]],cex=1.1)

res = two_d_plot_visualize_covariate(d$PC5[inds],
    d$PC6[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC5",ylab="PC6",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="topright",names(res[[1]]),
       fill = res[[1]],cex=1.1)

manual_clustering = d$PC1 < 0 & d$PC2 < 0.02
names(manual_clustering) = rownames(d)
table(manual_clustering,to_rem)
table(manual_clustering,d$Cohort)
table(manual_clustering,is_jap[names(manual_clustering)])

save(manual_clustering,file="manual_clustering.RData")

res = two_d_plot_visualize_covariate(d$PC1[inds],
  d$PC2[inds],manual_clustering[inds],manual_clustering[inds],
  main = "Inferred EUs",xlab="PC1",ylab="PC2",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="topright",names(res[[1]]),fill = res[[1]],cex=1.1,ncol = 1)

inds = 1:nrow(d)
cohort_name = c("Cooper","ELITE")
res = two_d_plot_visualize_covariate(d$PC13[inds],
    d$PC12[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC13",ylab="PC12",lwd=2,cex.axis=1.4,cex.lab=1.4,
    xlim = c(-0.1,0.1),ylim = c(-0.1,0.1))
legend(x="topleft",cohort_name[as.numeric(names(res[[1]]))],
       fill = res[[1]],cex=1.3)

res = two_d_plot_visualize_covariate(d$PC2[inds],
      d$PC3[inds],kmeans_res[inds],kmeans_res[inds],
      main = "PCA+Clustering",xlab="PC2",ylab="PC3",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="bottomright",names(res[[1]]),fill = res[[1]],cex=1.3)

res = two_d_plot_visualize_covariate(d$PC3[inds],
    d$PC2[inds],d$Cohort[inds],d$Cohort[inds],
    main = "Cooper and Elite",xlab="PC3",ylab="PC2",lwd=2,cex.axis=1.4,cex.lab=1.4)
legend(x="topleft",cohort_name[as.numeric(names(res[[1]]))],
       fill = res[[1]],cex=1.3)

res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],d$batch[inds],d$batch[inds],
    main = "Cooper and Elite",xlab="PC1",ylab="PC2",cex.axis=1.4,cex.lab=1.4)
legend(x="topleft",cohort_name[as.numeric(names(res[[1]]))],
       fill = res[[1]],cex=1.3)

curr_is_jap = is_jap[rownames(d)]
res = two_d_plot_visualize_covariate(d$PC1[inds],
    d$PC2[inds],curr_is_jap[inds],curr_is_jap[inds],
    main = "Is JA sample?",xlab="PC1",ylab="PC2")
legend(x="topleft",cohort_name[as.numeric(names(res[[1]]))],
       fill = res[[1]],cex=1.3)

# Check number of clusters in PCA plot
wss <- sapply(1:10,
              function(k){kmeans(pc_x, k, nstart=50,iter.max = 15 )$tot.withinss})
wss <- sapply(1:20,function(k){tot_wss_hluct(k,h,pc_x)})

plot(1:length(wss), wss,
     type="b", pch = 19, frame = FALSE,cex.lab=1.5,lwd=2,cex.axis=1.4,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

pc_ps = c()
for(j in 1:20){
  # pc_ps[j] = compute_pc_vs_discrete_variable_association_p(d[,paste("PC",j,sep="")],
  #                                                          d$Cohort)
  pc_ps[j] = compute_pc_vs_binary_variable_association_p(d[,paste("PC",j,sep="")],
                                                         d$Cohort)
}
p.adjust(pc_ps)


pc_rocs = c()
for(j in 1:40){
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="Cooper"
  p1 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="genepool"
  p2 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "Cooper" | d[,"Cohort"]=="genepool"
  p3 = compute_pc_vs_binary_variable_association_roc(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  pc_rocs = rbind(pc_rocs,c(p1,p2,p3))
}
pc_rocs > 0.75
pc_ps=c()
for(j in 1:40){
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="Cooper"
  p1 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "ELITE" | d[,"Cohort"]=="genepool"
  p2 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  curr_inds = d[,"Cohort"] == "Cooper" | d[,"Cohort"]=="genepool"
  p3 = compute_pc_vs_binary_variable_association_p(
    pc = pcax[curr_inds,paste("PC",j,sep="")],y = d[curr_inds,"Cohort"]
  )
  pc_ps = rbind(pc_ps,c(p1,p2,p3))
}
pc_ps < 1e-5
