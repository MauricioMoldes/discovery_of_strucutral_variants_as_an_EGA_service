###########################
## Coverage and Insert size 
###########################

cat ../data/$1/analysis/$1.depth_of_coverage | cut -d"=" -f2 | xargs  > ../data/$1/analysis/depth_coverage
cat ../data/$1/analysis/$1.stats |  grep "insert size average" |  cut -d":" -f2 | xargs  > ../data/$1/analysis/insert_size

echo "$1" > ../data/$1/analysis/sample_id

paste  ../data/$1/analysis/sample_id ../data/$1/analysis/depth_coverage ../data/$1/analysis/insert_size  > ../data/$1/analysis/join_depth_coverage_insert_size

cat ../data/sample_metrics.header ../data/$1/analysis/join_depth_coverage_insert_size > ../data/$1/analysis/sample_metrics

rm ../data/$1/analysis/join_depth_coverage_insert_size 
rm ../data/$1/analysis/depth_coverage
rm ../data/$1/analysis/insert_size
rm ../data/$1/analysis/$1.flagstat.statistics
rm ../data/$1/analysis/$1.stats
rm ../data/$1/analysis/$1.depth_of_coverage
rm ../data/$1/analysis/sample_id
