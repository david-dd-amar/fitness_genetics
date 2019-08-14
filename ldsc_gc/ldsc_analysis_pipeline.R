
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
  
}




# Prepare the input files
system(paste("mkdir ld_hub_input"))
for(file in allfiles){
  print(file)
  d = fread(file,data.table = F,stringsAsFactors = F)
  currbim = bim[d$SNP,]
  inds = is.element(d$SNP,set = names(rsids))
  d$SNP[inds] = rsids[d$SNP[inds]]
  ldhub_out = data.frame(
    snpid = d$SNP,
    A1 = d$A1,
    A2 = currbim$V6,
    N = d$NMISS,
    Zscore = d$STAT,
    Pvalue = d$P
  )
  colnames(ldhub_out) = c("SNP","A1","A2","N","Z","P")
  outfile = paste("ld_hub_input/",file,sep="")
  write.table(ldhub_out,file=outfile,row.names = F,col.names = T,sep="\t",quote = F)
  # system(paste("zip",paste(outfile,".zip",sep=""),outfile))
}

for(file in allfiles){
  outfile = paste("ld_hub_input/",file,sep="")
  system(paste("zip",paste(outfile,".zip",sep=""),outfile,"&"))
}

system("source activate ldsc")
system("python ../ldsc/ldsc.py --h2 cooper_treadmill_time_PCs4.assoc.linear --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cooper_treadmill")

# Some commands
# python ../ldsc/munge_sumstats.py --sumstats elite_vo2_ml_kg_min_PCs4.assoc.linear --merge-alleles w_hm3.snplist --out elite_vo2
# python ../ldsc/munge_sumstats.py --sumstats elite_vs_gp_gwas_res_PCs4.assoc.logistic --merge-alleles w_hm3.snplist --out elite_gp
# python ../ldsc/munge_sumstats.py --sumstats cooper_vs_gp_gwas_res_PCs4.assoc.logistic --merge-alleles w_hm3.snplist --out cooper_gp
# python ../ldsc/munge_sumstats.py --sumstats cooper_treadmill_time_PCs4.assoc.linear --merge-alleles w_hm3.snplist --out cooper_treadmill
# python ../ldsc/ldsc.py --rg elite_vo2.sumstats.gz,cooper_treadmill.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cooper_treadmill_vs_elite_vo2
# python ../ldsc/ldsc.py --rg elite_gp.sumstats.gz,cooper_gp.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cooper_vs_elite_gp









