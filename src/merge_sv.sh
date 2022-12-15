echo "../data/$1/lumpy_out/lumpy.vcf" >> ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss 
echo "../data/$1/manta_out/manta.vcf" >> ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss
echo "../data/$1/lumpy_out/lumpy.vcf" >> ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss 
echo "../data/$1/gridss_out/gridss.vcf" >> ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss

# by default sv-callers output (all.vcf) is the sum of the 4 different tools 
# here, we create a new merged vcf with the following rules : 

# Number supporting callers - 2,3 and 4
# SV Distance - 100
# Same strain - N  
# Min SV size - 50

../bin/SURVIVOR/Debug/./SURVIVOR merge ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss 100 2 1 1 1 50 ../data/$1/analysis/2.vcf
../bin/SURVIVOR/Debug/./SURVIVOR merge ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss 100 3 1 1 1 50 ../data/$1/analysis/3.vcf
../bin/SURVIVOR/Debug/./SURVIVOR merge ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss 100 4 1 1 1 50 ../data/$1/analysis/4.vcf

rm ../data/$1/analysis/sample_files_manta_delly_lumpy_gridss
