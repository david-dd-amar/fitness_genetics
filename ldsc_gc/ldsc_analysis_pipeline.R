
library(data.table,lib.loc = "~/R/packages")

# direct data
is_direct_geno = T
gwas_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/eu_gwas"
bimfile = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/merged_mega_data_autosomal.bim"
check_bim_map_to_rsids = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_with_genepool_/1000g/ID-merged_mega_data_autosomal-1000G.txt"

# imputed data
is_direct_geno = F
gwas_dir = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_reclustered_imp/eu_gwas_maf0.01/gwas_num_pcs_5/"
bimfile = NULL
check_bim_map_to_rsids = NULL

setwd(gwas_dir)
ldsc_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/ldsc/ldsc/"
ldsc_snp_list_path = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/ldsc/from_toturial/"
ldsc_snps = fread(paste(ldsc_snp_list_path,"w_hm3.snplist",sep=""),stringsAsFactors = F,data.table = F)

if(is_direct_geno){
  allfiles = list.files(".")
  allfiles = allfiles[grepl("linear$",allfiles)|grepl("logistic$",allfiles)]
  allfiles = allfiles[!grepl("fuma",allfiles)]
  
  # Set the required data to map snp ids to rsids
  bim = fread(bimfile,data.table = F,stringsAsFactors = F)
  rownames(bim) = bim[,2]
  bim_to_rsid = fread(check_bim_map_to_rsids,data.table = F,stringsAsFactors = F)
  rsids = sapply(bim_to_rsid[,2],function(x)strsplit(x,split=":")[[1]][1])
  names(rsids) = bim_to_rsid[,1]
  table(grepl("rs",rsids))
}
if(!is_direct_geno){
  allfiles = list.files(".")
  allfiles = allfiles[grepl("all",allfiles)]
  allfiles = allfiles[!grepl("fuma",allfiles)]
}

# Prepare the input files
system(paste("mkdir ld_hub_input"))
setwd(gwas_dir)
for(file in allfiles){
  print(file)
  d = fread(file,data.table = F,stringsAsFactors = F)
  if(is_direct_geno){
    currbim = bim[d$SNP,]
    inds = is.element(d$SNP,set = names(rsids))
    d$SNP[inds] = rsids[d$SNP[inds]]
    d$A2 = currbim$V6
  }
  zstat_ind = which(grepl("stat",names(d),ignore.case = T))[1]
  N_ind = which(grepl("nmiss|obs_ct",names(d),ignore.case = T))[1]
  p_ind = which(names(d)=="P")[1]
  rsid_ind = which(grepl("ID|SNP",names(d),ignore.case = F))[1]
  a2_ind = which(grepl("ALT|A1",names(d),ignore.case = F))[1]
  a1_ind = which(grepl("REF|A2",names(d),ignore.case = F))[1]
  ldhub_out = d[,c(rsid_ind,a1_ind,a2_ind,N_ind,zstat_ind,p_ind)]
  colnames(ldhub_out) = c("SNP","A1","A2","N","Z","P")
  rows_to_rem = is.na(as.numeric(as.character(ldhub_out$Z)))
  ldhub_out = ldhub_out[!rows_to_rem,]
  
  # Optional: reduce to the ldsc snp set
  ldhub_out = ldhub_out[is.element(ldhub_out$SNP,set=ldsc_snps$SNP),]
  print("ldsc input data size:")
  print(dim(ldhub_out))
  
  outfile = paste("ld_hub_input/",file,sep="")
  write.table(ldhub_out,file=outfile,row.names = F,col.names = T,sep="\t",quote = F)
  # system(paste("zip",paste(outfile,".zip",sep=""),outfile))
}

# for(file in allfiles){
#   outfile = paste("ld_hub_input/",file,sep="")
#   system(paste("zip",paste(outfile,".zip",sep=""),outfile,"&"))
# }

# make sure this is run outside R (?)
system("source activate ldsc")
setwd(paste(gwas_dir,"ld_hub_input/",sep=""))
for(file in allfiles){
  system(paste(
    "python", paste(ldsc_path,"munge_sumstats.py",sep=""),
    "--sumstats", paste(gwas_dir,"ld_hub_input/",file,sep=""),
    "--merge-alleles", paste(ldsc_snp_list_path,"w_hm3.snplist",sep=""),
    "--out",file
  ))
  system(paste(
    "python", paste(ldsc_path,"ldsc.py",sep=""),
    "--h2", paste(file,".sumstats.gz",sep=""),
    "--ref-ld-chr", paste(ldsc_snp_list_path,"eur_w_ld_chr/",sep=""),
    "--w-ld-chr", paste(ldsc_snp_list_path,"eur_w_ld_chr/",sep=""),
    "--out",paste(file,"_h2",sep="")
  ))
}

# Some commands
# python ../ldsc/munge_sumstats.py --sumstats elite_vo2_ml_kg_min_PCs4.assoc.linear --merge-alleles w_hm3.snplist --out elite_vo2
# python ../ldsc/munge_sumstats.py --sumstats elite_vs_gp_gwas_res_PCs4.assoc.logistic --merge-alleles w_hm3.snplist --out elite_gp
# python ../ldsc/munge_sumstats.py --sumstats cooper_vs_gp_gwas_res_PCs4.assoc.logistic --merge-alleles w_hm3.snplist --out cooper_gp
# python ../ldsc/munge_sumstats.py --sumstats cooper_treadmill_time_PCs4.assoc.linear --merge-alleles w_hm3.snplist --out cooper_treadmill
# python ../ldsc/ldsc.py --rg elite_vo2.sumstats.gz,cooper_treadmill.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cooper_treadmill_vs_elite_vo2
# python ../ldsc/ldsc.py --rg elite_gp.sumstats.gz,cooper_gp.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cooper_vs_elite_gp









