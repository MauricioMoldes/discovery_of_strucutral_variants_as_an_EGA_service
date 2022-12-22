################
## Number events chromossome
################

cat ../data/$1/analysis/all.decoy.size.vcf |  grep -v "#" |  cut -f1 > ../data/$1/analysis/chr
cat ../data/$1/analysis/all.decoy.size.vcf|  grep -v "#" |  cut -f8 | cut -d";" -f4 | cut -d "=" -f2  > ../data/$1/analysis/sv_type

paste ../data/$1/analysis/chr ../data/$1/analysis/sv_type |  sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' | grep -v "hs37d5" |  grep -v "GL*" | awk -v OFS="\t" '$1=$1'  > ../data/$1/analysis/tmp_sv_by_chr

cat ../data/number_events_chr.header  ../data/$1/analysis/tmp_sv_by_chr > ../data/$1/analysis/sv_by_chr

rm ../data/$1/analysis/chr
rm ../data/$1/analysis/sv_type
rm ../data/$1/analysis/tmp_sv_by_chr