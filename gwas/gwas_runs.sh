bfile=/oak/stanford/groups/euan/projects/fitness_genetics/analysis/mega_all_batches/merged_mega_data
pheno=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/master_phe_mega.phe
pcs="EU_PC1,EU_PC2,EU_PC3,EU_PC4,EU_PC5"
pop=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/eu_cluster1.phe
maf=0.05
excl=/oak/stanford/groups/euan/projects/fitness_genetics/pheno/rel_to_remove.phe
vexcl=/oak/stanford/groups/euan/projects/fitness_genetics/analysis/all_bad_snps.txt

ml load plink2

plink2 --bfile ${bfile} --chr 1-22 \
	--pheno ${pheno} --pheno-name Treadmill.time.sec \
	--covar ${pheno} --covar-name sex,${pcs} \
	--keep ${pop} --remove ${excl} \
	--glm skip hide-covar --maf ${maf} --adjust \
	--exclude ${vexcl} \
	--require-covar --require-pheno \
	--out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name VO2max..ml.kg.min. \
        --covar ${pheno} --covar-name sex,${pcs} \
        --keep ${pop}  --remove ${excl} \
	--exclude ${vexcl} \
        --glm skip hide-covar --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out elite

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name cooper_vs_elite \
        --covar ${pheno} --covar-name sex,${pcs} \
        --keep ${pop}  --remove ${excl} \
        --exclude ${vexcl} \
	--logistic skip hide-covar --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name cooper_vs_gp \
        --covar ${pheno} --covar-name sex,${pcs} \
        --keep ${pop}  --remove ${excl} \
        --logistic skip hide-covar --maf ${maf} --adjust \
        --exclude ${vexcl} \
	--require-covar --require-pheno \
        --out cooper

plink2 --bfile ${bfile} --chr 1-22 \
        --pheno ${pheno} --pheno-name elite_vs_gp \
        --covar ${pheno} --covar-name sex,${pcs} \
        --keep ${pop}  --remove ${excl} \
        --exclude ${vexcl} \
	--logistic skip hide-covar --maf ${maf} --adjust \
        --require-covar --require-pheno \
        --out elite
