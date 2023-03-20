
###############
## Number events per sample
################

cat ../data/$1/analysis/all.decoy.size.vcf | grep -v "#"  | cut -f8 | cut -d";" -f4 | cut -d "=" -f2 | sort | uniq -c |   sed -e 's/^[[:space:]]*//'  | awk -v OFS="\t" '$1=$1' | sed "s/^/$1\t/g" > ../data/$1/analysis/number_events_sv_no_header
cat ../data/number_events_sv.header ../data/$1/analysis/number_events_sv_no_header > ../data/$1/analysis/number_events_sv

rm ../data/$1/analysis/number_events_sv_no_header
