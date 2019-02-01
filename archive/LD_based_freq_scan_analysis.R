# In this input type we take a BIM and a result of check_bim from Wryner.
# The goal is to compare LD-based patterns to a reference. The assumption is that true deviations
# from the refenrece come in blocks.

freqplot_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/FreqPlot-merged_mega_data_autosomal-HRC.txt"
bim_file = "/oak/stanford/groups/euan/projects/fitness_genetics/analysis/no_recl_mega_separate_recalls/merged_mega_data_autosomal-hrc_updated.bim"

freq_info = read.table(freqplot_file,header = F,stringsAsFactors = F,row.names = 1)
bim_info = read.table(bim_file,header=F,stringsAsFactors = F,row.names = 2)

script_file = "/oak/stanford/groups/euan/projects/fitness_genetics/scripts/fitness_genetics/R/gwas_flow_helper_functions.R"
source(script_file)

d = bim_info[bim_info[,1]!="0",]
ord = order(d[,1],d[,3])
cor(1:length(ord),ord)> 0.9999

inds = intersect(rownames(bim_info),rownames(freq_info))
d = bim_info[inds,]
ord = order(d[,1],d[,3])
cor(1:length(ord),ord)> 0.9999

maf_diffs = freq_info[inds,3]
acf(maf_diffs,lag.max = 5,plot=F)
acf(freq_info[,2],lag.max = 5,plot=F)

table(maf_diffs>0.01)
x = as.numeric(maf_diffs>0.01)[-1]
g_dists = bim_info[inds[-1],3] - bim_info[inds[-length(inds)],3]
chisq.test(table(x,g_dists>10000))
ws = seq(1,length(x),by=10)
check_window<-function(i,by,x,num_1s = 2){
  is_1 = x[i]
  i1 = max(1,i-by)
  i2 = min(length(x),i+by)
  x = x[i1:i2]
  return(!is_1 || (is_1 && sum(x)>num_1s))
}
check_w_res = sapply(1:length(x),check_window,by=5,x=x)
table(check_w_res)
