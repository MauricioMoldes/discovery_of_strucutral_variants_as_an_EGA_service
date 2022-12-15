
cat ../data/$1/analysis/size_distribution_events |  tail -n +2 | sort | uniq -c |  sort -nr | sed 's/^ *//g' | sed 's/ /'$'\t''/' > ../data/$1/analysis/sv_size_discovered_no_header
cat ../data/size_count_distribution_events.header ../data/$1/analysis/sv_size_discovered_no_header >  ../data/$1/analysis/sv_size_discovered

rm ../data/$1/analysis/sv_size_discovered_no_header
