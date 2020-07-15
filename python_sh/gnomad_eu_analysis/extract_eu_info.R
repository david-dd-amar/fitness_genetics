args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
	print("usage:<gnomad vcf><out file>")
	q(save="no")
}

library(data.table)
vcf = args[1]
d = fread(vcf,stringsAsFactors=F,data.table=F,skip = "CHROM")
info =sapply(d[,8],function(x)strsplit(x,split=";")[[1]])

get_vals<-function(x,reg="fin",vnames=c("AC","AN","AF","nhomalt")){
	vnames = paste0(vnames,"_",reg)
	reg = paste0("_",reg,"=")
	inds = grepl(reg,x)
	v = x[inds]
	l = strsplit(v,split="=")
	newv = sapply(l,function(x)x[2])
	names(newv) = sapply(l,function(x)x[1])
	newv2 = rep(NA,length(vnames))
	names(newv2) = vnames
	intrs = intersect(names(newv),names(newv2))
	newv2[intrs] = newv[intrs]
	return(newv2)
}

finnish_info = sapply(info,get_vals,reg="fin")
non_finnish_eu_info = sapply(info,get_vals,reg="nfe")

colnames(finnish_info) = NULL
finnish_info = t(finnish_info)
colnames(non_finnish_eu_info) = NULL
non_finnish_eu_info = t(non_finnish_eu_info)
info=NULL
gc()

info = cbind(finnish_info,non_finnish_eu_info)
d = cbind(d[,-8],info)
newid = paste(gsub("chr","",d[,1]),d[,2],d[,4],d[,5],sep=":")
d = cbind(newid,d)
fwrite(d,file=args[2],row.names=F,col.names=T,sep="\t",quote=F)
