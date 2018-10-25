# Here we interpret the results of PCAngsd

dir = "/Users/David/Desktop/elite/sept2018_prepro_res/pcangsd/"
prefix = "cooper"

setwd(dir)

chisq_scores = as.matrix(read.table(paste(prefix,".selection",sep="")))
rownames(chisq_scores) = read.table(paste(prefix,".sites",sep=""),stringsAsFactors = F)[,1]
chisq_pvals = pchisq(chisq_scores,1,lower.tail = F)

apply(chisq_pvals<1e-10,2,sum)
qqplot(-log(rnorm(1000)),-log(chisq_pvals[sample(1:nrow(chisq_pvals))[1:10000],3]));abline(0,1)

x = cbind(rownames(chisq_pvals),chisq_pvals[,2])
colnames(x) = c("rsID","P-value")
write.table(x,file=paste(prefix,"_pc2_fuma.txt",sep=""),row.names = F,col.names = T,quote = F,sep=" ")
