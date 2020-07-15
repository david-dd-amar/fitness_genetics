for i in {1..22}
do
Rscript compare_freq_with_gnomad.R gnomad_intersect/chr${i}.eu_info.tsv eu_pgen/chr${i}.afreq comp_freq_with_gnomad.chr${i}.tsv > comp_freq_with_gnomad.chr${i}.log  &
done

i=1
less comp_freq_with_gnomad.chr${i}.tsv > comp_freq_with_gnomad.all.tsv
for i in {2..22}
do
less comp_freq_with_gnomad.chr${i}.tsv >> comp_freq_with_gnomad.all.tsv
done
