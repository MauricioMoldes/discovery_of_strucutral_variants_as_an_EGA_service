###############
## index vcf 
###############

bgzip -c ../data/$1/lumpy_out/lumpy.vcf > ../data/$1/lumpy_out/lumpy.vcf.gz
tabix -p vcf ../data/$1/lumpy_out/lumpy.vcf.gz

bgzip -c ../data/$1/manta_out/manta.vcf > ../data/$1/manta_out/manta.vcf.gz
tabix -p vcf ../data/$1/manta_out/manta.vcf.gz

bgzip -c ../data/$1/delly_out/delly.vcf > ../data/$1/delly_out/delly.vcf.gz
tabix -p vcf ../data/$1/delly_out/delly.vcf.gz

bgzip -c ../data/$1/gridss_out/lumpy.vcf > ../data/$1/gridss_out/gridss.vcf.gz
tabix -p vcf ../data/$1/gridss_out/gridss.vcf.gz