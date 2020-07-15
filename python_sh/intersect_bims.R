library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
	print("usage:<bim1><bim2><outpath>")
	q(save="no")
}
print(paste("intersecting two bim files:"))
bim1 = args[1]
bim2 = args[2]
# Read the bim files
b1 = fread(bim1,stringsAsFactors=F,data.table=F)
print(paste("File 1:",bim1))
print(head(b1))
b2 = fread(bim2,stringsAsFactors=F,data.table=F)
print(paste("File 2:",bim2))
print(head(b2))
# change the snp ids to be in the same format
# we assume that these bed files were created 
# using the reference genome as the major allele
b1_original_ids = b1[,2]
b2_original_ids = b2[,2]
b1[,2] = paste(b1[,1],b1[,4],b1[,6],b1[,5],sep=":")
b2[,2] = paste(b2[,1],b2[,4],b2[,6],b2[,5],sep=":")
# print the new bim files
bim1 = gsub(".bim",".new.bim",bim1)
bim2 = gsub(".bim",".new.bim",bim2)
fwrite(b1,file=bim1,col.names=F,row.names=F,sep="\t",quote=F)
fwrite(b2,file=bim2,col.names=F,row.names=F,sep="\t",quote=F)

print("id intersect:")
ids = intersect(b1[,2],b2[,2])
print(length(ids))
fwrite(t(t(ids)),file=paste0(args[3],"intersect_variants.txt"),quote=F,row.names=F,col.names=F)
ref_info = b1[b1[,2] %in% ids,c(2,6)]
fwrite(ref_info,file=paste0(args[3],"intersect_variants_ref.txt"),quote=F,row.names=F,col.names=F,sep="\t")

print("checking other intersections for QC")
print("Intersection of original ids:")
print(length(intersect(b2_original_ids,b1_original_ids)))
print("Intersect of alternative ids - chr:pos:col5:col6")
b1_ids_alt = paste(b1[,1],b1[,4],b1[,5],b1[,6],sep=":")
print(length(intersect(b2[,2],b1_ids_alt)))
