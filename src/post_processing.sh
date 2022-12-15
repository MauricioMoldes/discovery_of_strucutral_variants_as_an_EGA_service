##########################
## Post-Processing steps
##########################

##########################
## Filter PASS flag
##########################

#bcftools view -f PASS ../data/$1/lumpy_out/lumpy.vcf > ../data/$1/analysis/lumpy.pass.vcf
#bcftools view -f PASS ../data/$1/manta_out/manta.vcf > ../data/$1/analysis/manta.pass.vcf
#bcftools view -f PASS ../data/$1/lumpy_out/delly.vcf > ../data/$1/analysis/delly.pass.vcf
#bcftools view -f PASS ../data/$1/gridss_out/gridss.vcf > ../data/$1/analysis/gridss.passvcf

##########################
## Filter decoy sequences variants
##########################

#bcftools view ../data/$1/analysis/lumpy.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/lumpy.pass.decoy.vcf
#bcftools view ../data/$1/analysis/manta.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/manta.pass.decoy.vcf
#bcftools view ../data/$1/analysis/delly.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/delly.pass.decoy.vcf
#bcftools view ../data/$1/analysis/gridss.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/gridss.pass.decoy.vcf

##########################
## Filter Genotype unknown or ref-ref
##########################



##########################
## Filter Read depth 10 < SV < 1000
##########################



##########################
## Filter size 50 < SV < 10000
##########################

