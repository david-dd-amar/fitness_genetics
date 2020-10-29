bfile=/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_all_batches/merged_mega_data
pheno=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/master_phe_mega.phe
pcs="EU_PC1,EU_PC2,EU_PC3,EU_PC4,EU_PC5"
maf=0.01
excl1=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/rel_to_remove.phe
excl2=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_pca_outliers.phe
vexcl=/oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt
cols_add=cols=+ax,+a1countcc,+gcountcc,+a1freq,+a1freqcc,+a1count,+totallele

excl="tmp_all_samps_to_rem.txt"
less ${excl1} > ${excl}
less ${excl2} >> ${excl}

ml load plink2

plink2 --bfile ${bfile} --chr 1-22 \
	--pheno ${pheno} --pheno-name Treadmill.time.sec \
	--covar ${pheno} --covar-name sex,${pcs} \
	--remove ${excl} \
	--glm skip hide-covar ${cols_add} --maf ${maf} --adjust \
	--exclude ${vexcl} \
	--require-covar --require-pheno \
	--out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name VO2max..ml.kg.min. \
        --covar ${pheno} --covar-name sex,${pcs} \
        --remove ${excl} \
	--exclude ${vexcl} \
        --glm skip hide-covar ${cols_add} --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out elite

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name cooper_vs_elite \
        --covar ${pheno} --covar-name sex,${pcs} \
        --remove ${excl} \
        --exclude ${vexcl} \
	--logistic skip hide-covar ${cols_add} --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name cooper_vs_gp \
        --covar ${pheno} --covar-name sex,${pcs} \
        --remove ${excl} \
        --logistic skip hide-covar ${cols_add} --maf ${maf} --adjust \
        --exclude ${vexcl} \
	--require-covar --require-pheno \
        --out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name elite_vs_gp \
        --covar ${pheno} --covar-name sex,${pcs} \
        --remove ${excl} \
        --exclude ${vexcl} \
	--logistic skip hide-covar ${cols_add} --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out elite
