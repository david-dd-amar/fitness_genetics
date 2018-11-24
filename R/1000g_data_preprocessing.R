wd = "/oak/stanford/groups/euan/projects/fitness_genetics/1000g/"
setwd(wd)

script_file = "~/repos/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

curr_gz_files = list.files()
curr_gz_files = curr_gz_files[grep("vcf.gz$",curr_gz_files)]
curr_names = sapply(curr_gz_files,function(x)strsplit(x,split='\\.')[[1]][2])
names(curr_gz_files) = curr_names

# Preprocessing
# Flipscan each chromosome
for(nn in curr_names){
  curr_cmd = paste("plink2 --vcf",curr_gz_files[nn],
                   "--maf 0.001",
                   "--threads 4",
                   "--make-bed --out",nn
  )
  run_plink_command(curr_cmd,wd,nn,
                    get_sh_prefix_one_node_specify_cpu_and_mem,Ncpu=4,mem_size=32000,
                    plink_pkg = "plink/2.0a1")
}
wait_for_job(waittime = 120)

# Correct the fam files
for(chr in 1:22){
  fam = read.table(paste("chr",chr,".fam",sep=""),stringsAsFactors = F)
  fam[,1] = fam[,2]
  write.table(fam,paste("chr",chr,".fam",sep=""),sep="\t",row.names = F,col.names = F,
              quote=F)
}
