############################
### ANNOTATIONS
############################

java -Xmx4g -jar ../bin/snpEff/snpEff.jar -v GRCh37.75 ../data/$1/analysis/all.vcf > ../data/$1/analysis/all.annot.vcf

cat ../data/$1/analysis/all.annot.vcf | grep -v "#" | cut -f8 | cut -d";" -f11 | cut -d"|" -f2 | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//'  | awk -v OFS="\t" '$1=$1' | awk -v FS='\t' -v OFS='\t' '{$2=$2!=""?$2:"NA"}1'  >   ../data/$1/analysis/annotation_counts_no_header 
cat ../data/annotations.header ../data/$1/analysis/annotation_counts_no_header > ../data/$1/analysis/annotation_counts

rm ../data/$1/analysis/annotation_counts_no_header
