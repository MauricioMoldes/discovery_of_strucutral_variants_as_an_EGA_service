##################
## Create sample folder
##################

mkdir -p ../data/$1 ../data/$1/analysis ../data/$1/delly_out ../data/$1/manta_out ../data/$1/gridss_out ../data/$1/lumpy_out

################################
## Pull vcf data from AE01
#################################

scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/$1/all.vcf  ../data/$1/analysis/all.vcf
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/$1/lumpy_out/lumpy.vcf  ../data/$1/lumpy_out/lumpy.vcf 
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/$1/manta_out/manta.vcf  ../data/$1/manta_out/manta.vcf 
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/$1/delly_out/delly.vcf  ../data/$1/delly_out/delly.vcf 
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/$1/gridss_out/gridss.vcf ../data/$1/gridss_out/gridss.vcf

##############################
## Pull statistics from AE01
##############################

scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/statistics/$1.flagstat.statistics  ../data/$1/analysis/$1.flagstat.statistics
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/statistics/$1.stats  ../data/$1/analysis/$1.stats
scp -oProxyCommand="ssh -W %h:%p mmoldes@gw-in.ega.crg.eu" devers@apps-executor-01u:/vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/statistics/$1.depth_of_coverage  ../data/$1/analysis/$1.depth_of_coverage

