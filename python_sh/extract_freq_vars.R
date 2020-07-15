library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=4){
        print("usage:<freq1><freq2><thr1><thr2>")
        q(save="no")
}

fr1 = fread(args[1],stringsAsFactors=F,data.table=F)
fr2 = fread(args[2],stringsAsFactors=F,data.table=F)
thr1 = as.numeric(args[3])
thr2 = as.numeric(args[4])
mafs1 = pmin(fr1[["ALT_FREQS"]],1-fr1[["ALT_FREQS"]])
mafs2 = pmin(fr2[["ALT_FREQS"]],1-fr2[["ALT_FREQS"]])

fr1 = fr1[mafs1 >= thr1,]
fr2 = fr2[mafs2 >= thr2,]

vars = intersect(fr1[["ID"]],fr2[["ID"]])
write.table(t(t(vars)),quote=F,row.names=F,col.names=F)
