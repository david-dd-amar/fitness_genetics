args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
        print("usage:<phe file with group column><king table><out file name>")
        q("no")
}

d1 = read.table(args[1],header=T,stringsAsFactors=F)
d2 = read.table(args[2],stringsAsFactors=F)
head(d1)
head(d2)
rownames(d1) = as.character(d1[,1])

allsamps = unique(union(d2[,1],d2[,2]))
print(table(d1[allsamps,2]))

samps = c()
for(i in 1:nrow(d2)){
	s1 = d2[i,1]
	s2 = d2[i,2]
	c1 = d1[s1,2]
	c2 = d1[s2,2]
	if(c1 == c2){next}
	print(paste(s1,c1,s2,c2,sep=","))
	#if(c1 == c2){next}
	samps = unique(union(samps,c(s1,s2)))
}
print(paste("number of identified duplicates:",length(samps)))

write.table(t(t(samps)),file=args[3],quote=F,row.names=F,col.names=F)

