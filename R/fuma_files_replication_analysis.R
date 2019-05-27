setwd("/Users/David/Desktop/elite/april_2019/")

p_discovery = 5e-08
p_rep = 0.001

library(data.table)
get_replication_results<-function(discovery_path,replication_path,p_discovery){
  discovery = fread(discovery_path,sep=" ",data.table = F)
  discovery = discovery[!is.na(discovery[,3]),]
  discovery = discovery[discovery$`P-value`<= p_discovery,]
  discovery_set = paste(discovery[,1],discovery[,2])
  replication = fread(replication_path,sep=" ",data.table = F)
  rep_set = paste(replication[,1],replication[,2])
  discovery_set_in_rep = discovery_set[is.element(discovery_set,set=rep_set)]
  discovery_set_in_rep_ps = replication[is.element(rep_set,set=discovery_set_in_rep),]
  return(list(discovery_set_in_rep_ps=discovery_set_in_rep_ps,discovery_set=discovery_set))
}

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./special_analysis_vs_ukbb/fuma_cooper_gwas_res_all_pcs10.assoc"
rep_res = get_replication_results(discovery_path,replication_path,p_discovery)

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./special_analysis_vs_ukbb/fuma_cooper_gwas_res_all_pcs10.assoc"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./gp_vs_cooper_imp/fuma_cooper_gwas_res_all_pcs1.assoc"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./ukbb_activity/fuma_physical_activity.physical_activity.glm.linear"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)

discovery_path = "./special_analysis_vs_ukbb/fuma_cooper_gwas_res_10pcs_cleaned_1e6_0.001.assoc"
replication_path = "./ukbb_activity/fuma_physical_activity.physical_activity.glm.linear"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./ukbb_fitness/ehr_discovery_raw_filteredby_0.01.txt"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)

discovery_path = "./special_analysis_vs_ukbb/fuma_cooper_gwas_res_10pcs_cleaned_1e6_0.001.assoc"
replication_path = "./ukbb_fitness/ehr_discovery_raw.txt"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)
hist(rep_res$discovery_set_in_rep_ps[,3],breaks=100)

discovery_path = "./gp_vs_cooper/fuma_cooper_vs_gp_gwas_res.assoc.logistic"
replication_path = "./ukbb_fitness/ehr_discovery_raw.txt"
rep_res = get_replication_results(discovery_path,replication_path,5e-06)
hist(rep_res$discovery_set_in_rep_ps[,3])









