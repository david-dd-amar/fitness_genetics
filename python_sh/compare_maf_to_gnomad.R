args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=4){
	print("usage:<path with gnomad afs: tables, second column is variant id, there is an AF_nfe column><ref afreq file><max_abs_diff><max log OR>")
	q("no")
}
max_diff = as.numeric(args[3])
max_log_or = as.numeric(args[4])
chrs = 1:22
to_rem = c()
library(data.table)
ref_freqs = fread(args[2],data.table=F,stringsAsFactors=F)
rownames(ref_freqs) = ref_freqs[,2]
#head(ref_freqs)
for(chr in chrs){
	currf = paste0(args[1],".chr",chr,".tsv")
	currd = fread(currf,data.table=F,stringsAsFactors=F)
	#print(head(currd))
	#break
	rownames(currd) = currd[,2]
	shared = intersect(rownames(currd),rownames(ref_freqs))
	shared_mafs1 = ref_freqs[shared,"ALT_FREQS"]
	shared_mafs2 = currd[shared,"AF_nfe"]
	absdiffs = abs(shared_mafs1-shared_mafs2)
	abs_logOR = abs(log(shared_mafs1)-log(shared_mafs2))
	to_rem = c(to_rem,shared[absdiffs>max_diff | abs_logOR > max_log_or])
	print(paste("finished analyzing chr",chr,"total to remove:",length(to_rem)))
}

write.table(t(t(to_rem)),quote=F,row.names=F,col.names=F)
