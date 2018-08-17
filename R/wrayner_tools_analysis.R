# This script uses the tools from http://www.well.ox.ac.uk/~wrayner/tools/index.html on our datasets
# The goal is to do a rigorous QC and examination of the data we generated

# First, download and install the tools on the cluster
# wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# unzip HRC-1000G-check-bim-v4.2.9-NoReadKey.zip
# mkdir check_bim
# mv HRC* check_bim/
# mv LICENSE.txt check_bim/
# cd check_bim
# wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
# gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

bim_paths = list(
  "our_original" = "/oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering",
  "with_ukbb_simple_merge" = "",
  "with_ukbb_qctools" = "/oak/stanford/groups/euan/projects/ukbb/qctool/merged_bed_final_for_gwas"
)

check_bim_script = "/oak/stanford/groups/euan/projects/fitness_genetics/wrayner_tools/check_bim/HRC-1000G-check-bim-NoReadKey.pl"

# example for cmd:
# perl HRC-1000G-check-bim-NoReadKey.pl -b /oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.bim -f /oak/stanford/groups/euan/projects/fitness_genetics/illu_processed_plink_data/no_reclustering/PLINK_050618_0953/recall_may_2018_without_reclustering.frq -hrc -p EUR -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab