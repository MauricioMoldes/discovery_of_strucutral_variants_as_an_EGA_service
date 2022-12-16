##########################
## Post-Processing steps
##########################

##########################
## Filter PASS flag
##########################

bcftools view -f PASS ../data/$1/lumpy_out/lumpy.vcf.gz > ../data/$1/analysis/lumpy.pass.vcf.gz
bcftools view -f PASS ../data/$1/manta_out/manta.vcf.gz > ../data/$1/analysis/manta.pass.vcf.gz
bcftools view -f PASS ../data/$1/lumpy_out/delly.vcf.gz > ../data/$1/analysis/delly.pass.vcf.gz
bcftools view -f PASS ../data/$1/gridss_out/gridss.vcf.gz > ../data/$1/analysis/gridss.pass.vcf.gz

##########################
## Filter decoy sequences variants
##########################

bcftools view ../data/$1/analysis/lumpy.pass.vcf.gz  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/lumpy.pass.decoy.vcf.gz
bcftools view ../data/$1/analysis/manta.pass.vcf.gz | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/manta.pass.decoy.vcf.gz
bcftools view ../data/$1/analysis/delly.pass.vcf.gz  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/delly.pass.decoy.vcf.gz
bcftools view ../data/$1/analysis/gridss.pass.vcf.gz  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/gridss.pass.decoy.vcf.gz

##########################
## Filter Genotype unknown or ref-ref
##########################



##########################
## Filter Read depth 10 < SV < 1000
##########################



##########################
## Filter size 50 < SV < 10000
##########################



