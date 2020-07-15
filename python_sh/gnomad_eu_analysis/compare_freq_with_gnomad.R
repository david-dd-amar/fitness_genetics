args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
	print("usage:<gnomad file><freq file><out file>")
	q(save="no")
}

library(data.table)
#gnomad_file = "/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_filtered_vcfs/gnomad_intersect/chr21.eu_info.tsv"
#freq_file = "/oak/stanford/groups/euan/projects/elite/ukbb_exome/elite_filtered_vcfs/eu_pgen/chr21.afreq"

gnomad_file = args[1]
freq_file = args[2]
d1 = fread(gnomad_file,stringsAsFactors=F,data.table=F,skip = "CHROM")
d2 = fread(freq_file,stringsAsFactors=F,data.table=F,skip = "CHROM")
d2[["ID"]] = gsub("chr","",d2[["ID"]])

dups = names(which(table(d2[["ID"]])>1))
if(length(dups)>0){
	print(paste("removing",length(dups),"duplicated ids in freq file"))
	d2 = d2[!(d2[["ID"]] %in% dups),]
}

rownames(d2) = d2[["ID"]]
rownames(d1) = d1[,1]

shared = intersect(d2[["ID"]],d1[,1])
d1 = d1[shared,]
d2 = d2[shared,]

elite_ac = round(d2[,"ALT_FREQS"] * d2[,"OBS_CT"])
elite_non_ac = d2[,"OBS_CT"] - elite_ac

simple_prop_test<-function(x){
	ac = x[1]
	non_ac = x[2]
	p = x[3]
	if(is.na(ac)){return(NA)}
	if(is.na(non_ac)){return(NA)}
	if(is.na(p)){return(NA)}
	return(prop.test(as.table(c(ac,non_ac)),p=p)$p.value)
}

m = cbind(elite_ac,elite_non_ac,d1[,"AF_fin"])
m[,3] = pmax(m[,3],0.001)
m[,3] = pmin(m[,3],1-0.001)
finnish_ps = apply(m,1,simple_prop_test)

m = cbind(elite_ac,elite_non_ac,d1[,"AF_nfe"])
m[,3] = pmax(m[,3],0.001)
m[,3] = pmin(m[,3],1-0.001)
eu_ps = apply(m,1,simple_prop_test)

out = cbind(d2,d1[,"AF_fin"],finnish_ps,d1[,"AF_nfe"],eu_ps)
colnames(out)[ncol(d2)+1] = "AF_fin"
colnames(out)[ncol(d2)+2] = "P_fin"
colnames(out)[ncol(d2)+3] = "AF_nfe"
colnames(out)[ncol(d2)+4] = "P_nfe"
fwrite(out,file=args[3],row.names=F,sep="\t",quote=F,col.names=T)


