##########################
## Post-Processing steps
##########################

##########################
## Filter PASS flag
##########################

## already done by sv-callers survivor filter "pass"

##########################
## Filter decoy sequences variants
##########################

cat ../data/$1/analysis/all.vcf |  grep "#" > ../data/$1/analysis/all.header 
cat ../data/$1/analysis/all.vcf  | grep -v "#" |  grep -v "hs37d5" | grep  -v "GL000" > ../data/$1/analysis/all.decoy.variants
cat ../data/$1/analysis/all.header ../data/$1/analysis/all.decoy.variants > ../data/$1/analysis/all.decoy.vcf

rm ../data/$1/analysis/all.header
rm ../data/$1/analysis/all.decoy.variants

##########################
## Filter Genotype unknown or ref-ref
##########################

#TODO

##########################
## Filter Read depth 10 < SV < 1000
##########################

# TODO

##########################
## Filter size 50 bases < SV < 10.000 bases
##########################

cat ../data/$1/analysis/all.decoy.vcf | grep "#" > ../data/$1/analysis/all.header 
cat ../data/$1/analysis/all.decoy.vcf | grep -v "#" > ../data/$1/analysis/all.decoy.variants


input="../data/$1/analysis/all.decoy.variants"
while IFS= read -r line
do
   
  sv_type=$(echo "$line" | cut -f8 | cut -d";" -f3,4 | sed s/SVLEN=//g |  sed s/SVTYPE=//g |  sed s/-//g | cut -d";" -f2)
  sv_size=$(echo "$line" | cut -f8 | cut -d";" -f3,4 | sed s/SVLEN=//g |  sed s/SVTYPE=//g |  sed s/-//g | cut -d";" -f1)

  if [[ $sv_type != "TRA" ]]; then
  	if  (($sv_size >= 50 && $sv_size <= 10000)); then
 			echo "$line" >> ../data/$1/analysis/all.decoy.bases.vcf
 	fi
  else
  	echo "$line" >> ../data/$1/analysis/all.decoy.bases.vcf
  fi
done < "$input"

cat ../data/$1/analysis/all.header ../data/$1/analysis/all.decoy.bases.vcf > ../data/$1/analysis/all.decoy.size.vcf

rm ../data/$1/analysis/all.header
rm ../data/$1/analysis/all.decoy.variants
rm ../data/$1/analysis/all.decoy.bases.vcf


