# In this script we take two bed files and merge them as follows:
# Input is two bed paths and array annotation 
# 1. Pre-analysis: correct the strand issue, this is specific for the platform
# 2. We use William Rayner's bim check tool to align the data with HRC
#    This step may remove some SNPs
# 3. We use plink to transform the files to bgen
# 4. We merge the bgens using qctool
# 5. We convert the bgen to bed - the output
# 6. (Optional): given a case-control phenotype we check for issues using PLINK's flip scan

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


bfile1 = ""
annotation1 = ""
bfile2 = ""
annotation2 = ""
outf = ""

bimcheck_path = "/home/users/davidama/apps/check_bim/HRC-1000G-check-bim-NoReadKey.pl"
hrc_panel_info = "/home/users/davidama/apps/check_bim/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
qctool_path = "/home/users/davidama/apps/qctool_v2/build/release/qctool_v2.0.1"

bim_paths = list(
  "our_original" = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering",
  "with_ukbb_simple_merge" = "",
  "with_ukbb_qctools" = "/oak/stanford/groups/euan/projects/ukbb/qctool/merged_bed_final_for_gwas"
)

# example for cmd:
# perl HRC-1000G-check-bim-NoReadKey.pl -b /oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.bim -f /oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.frq -hrc -p EUR -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab