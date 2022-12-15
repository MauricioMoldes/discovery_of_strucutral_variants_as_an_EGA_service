cat ../data/$1/delly_out/delly.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_delly
cat ../data/$1/lumpy_out/lumpy.vcf  | grep -v "#" |  wc -l > ../data/$1/analysis/number_events_lumpy
cat ../data/$1/manta_out/manta.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_manta
cat ../data/$1/gridss_out/gridss.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_gridss
