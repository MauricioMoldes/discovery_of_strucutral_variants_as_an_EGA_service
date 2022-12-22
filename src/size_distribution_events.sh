cat ../data/$1/analysis/all.decoy.size.vcf  | grep -v "#"  | cut -f8 | cut -d";" -f3,4 | sed s/SVLEN=//g |  sed s/SVTYPE=//g | sed s/-//g | sed 's/;/\t/g' | sed "s/^/$1\t/g" > ../data/$1/analysis/size_distribution_events_no_header
cat ../data/size_distribution_events.header ../data/$1/analysis/size_distribution_events_no_header > ../data/$1/analysis/size_distribution_events

rm ../data/$1/analysis/size_distribution_events_no_header
