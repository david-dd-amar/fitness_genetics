# module load python/3.6.1 

import os
# os.system("module load py-pandas/0.23.0_py27")
# os.system("module load py-numpy/1.14.3_py27")
# os.system("module load py-scikit-image/0.15.0_py27")
import re
import gzip
import numpy as np
import pandas as pd

# input to the script
gz1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_cnv_data/penncnv_table_batch1.txt.gz"
gz2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_cnv_data/penncnv_table_other_batches.txt.gz"
snpfile1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/SNP_Table.txt"
snpfile2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/SNP_Table.txt"
sample_file1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/MEGA_Consortium_recall/Samples_Table.txt"
sample_file2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/non_MEGA_Cons_recall/Samples_Table.txt"
phe_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_with_genepool/integrated_sample_metadata_and_covariates.phe"

# output path
out1 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_cnv_data/penncnv_table_batch1_filtered/"
# os.system("mkdir "+out1)
out2 = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_cnv_data/penncnv_table_other_batches_filtered/"
# os.system("mkdir "+out2)

phe_data = pd.read_table(phe_file,sep=" ")
ids_in_phe = phe_data.IID.to_dict().values()

samples2 = pd.read_table(sample_file2)
samples2_set = samples2["Sample ID"].tolist()

snps1 = pd.read_table(snpfile1)
print("read in the first snp qc file, size:")
print(snps1.shape)
snps1 = snps1[snps1["Call Freq"]>0.95]
print("remaining snps after filtering:")
print(snps1.shape)
snps1_set = snps1["Name"].tolist()
snps1_set = set(snps1_set)

snps2 = pd.read_table(snpfile2)
print("read in the second snp qc file, size:")
print(snps2.shape)
snps2 = snps2[snps2["Call Freq"]>0.95]
print("remaining snps after filtering:")
print(snps2.shape)
snps2_set = snps2["Name"].tolist()
snps2_set = set(snps2_set)

shared_snps = snps1_set.intersection(snps2_set)
print("number of shared snps between the two MEGA platforms:")
print(len(shared_snps))
shared_snps_dict = {}
for snp in shared_snps:
    shared_snps_dict[snp] = True


# read the first line of the gz1 file
with gzip.open(gz1,'rt') as f:
    for line in f:
        gz1_header_arr = line.split("\t")
        break
    f.close()

# read the first line of the gz2 file
with gzip.open(gz2,'rt') as f:
    for line in f:
        gz2_header_arr = line.split("\t")
        break
    f.close()

def check_regexes(rs,x):
    for regex in rs:
        if re.search(regex,x):
            return True
    return False

# tests
# check_regexes(["aba","abba"],"dd")
# check_regexes(["aba","abba"],"nnn abba 888")

# Define the columns to ignore from the first batch:
# 1. samples that are not in the phe file, AND
# 2. samples that also appear in the gz2 file
gz1_inds_to_include = [1,3,4]
for i in range(len(gz1_header_arr)):
    currval = gz1_header_arr[i]
    if not check_regexes(ids_in_phe,currval):
        continue
    if check_regexes(samples2_set,currval):
        continue
    # print("keeping: " + gz1_header_arr[i])
    gz1_inds_to_include.append(i)

# Define the columns to ignore from the second batch:
# 1. samples that are not in the phe file, AND
# 2. There is no need to keep the Name, Chr, and position (we combine the two datasets)
gz2_inds_to_include = [1,3,4]
for i in range(len(gz2_header_arr)):
    currval = gz2_header_arr[i]
    if not check_regexes(ids_in_phe,currval):
        continue
    # print("keeping: " + gz1_header_arr[i])
    gz2_inds_to_include.append(i)

# Covered sample sizes
print("Finished filtering samples out, remaining number of columns (by file)")
print(len(gz1_inds_to_include)/3)
print(len(gz2_inds_to_include)/3)

all_chrs = snps1["Chr"].astype('str').unique()

o1s = {}
for chr in all_chrs:
    o1s[chr] = open(out1+chr+".txt","w")

with gzip.open(gz1,'rt') as f:
    for line in f:
        arr = line.split("\t")
        currChr = arr[3]
        # print(arr[1]+" " + arr[3]+" "+arr[4])
        if arr[1]!="Name" and not shared_snps_dict.has_key(arr[1]):
            continue
        arr_reduced = []
        for i in gz1_inds_to_include:
            arr_reduced.append(arr[i])
        s = "\t".join(arr_reduced)
        if arr[1]=="Name":
            for o1 in o1s.values():
                o1.write(s+"\n")
        else:
            o1s[currChr].write(s+"\n")
    f.close()
    for o1 in o1s.values():o1.close()

o2s = {}
for chr in all_chrs:
    o2s[chr] = open(out2+chr+".txt","w")

with gzip.open(gz2,'rt') as f:
    for line in f:
        arr = line.split("\t")
        currChr = arr[3]
        # print(arr[1]+" " + arr[3]+" "+arr[4])
        if arr[1]!="Name" and not shared_snps_dict.has_key(arr[1]):
            continue
        arr_reduced = []
        for i in gz2_inds_to_include:
            arr_reduced.append(arr[i])
        s = "\t".join(arr_reduced)
        if arr[1]=="Name":
            for o2 in o2s.values():
                o2.write(s+"\n")
        else:
            o2s[currChr].write(s+"\n")
    f.close()
    for o2 in o2s.values():o2.close()

print("completed analyzing the two files")




