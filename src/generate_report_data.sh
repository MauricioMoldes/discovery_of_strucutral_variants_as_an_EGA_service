#!/usr/bin/env bash

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

##########################
## Post-Processing steps
##########################

##########################
## Filter PASS flag
##########################

#bcftools view -f PASS ../data/$1/lumpy_out/lumpy.vcf > ../data/$1/analysis/lumpy.pass.vcf
#bcftools view -f PASS ../data/$1/manta_out/manta.vcf > ../data/$1/analysis/manta.pass.vcf
#bcftools view -f PASS ../data/$1/lumpy_out/delly.vcf > ../data/$1/analysis/delly.pass.vcf
#bcftools view -f PASS ../data/$1/gridss_out/gridss.vcf > ../data/$1/analysis/gridss.passvcf

##########################
## Filter decoy sequences variants
##########################

#bcftools view ../data/$1/analysis/lumpy.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/lumpy.pass.decoy.vcf
#bcftools view ../data/$1/analysis/manta.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/manta.pass.decoy.vcf
#bcftools view ../data/$1/analysis/delly.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/delly.pass.decoy.vcf
#bcftools view ../data/$1/analysis/gridss.pass.vcf  | grep -v "#" | grep -v "hs37d5" > ../data/$1/analysis/gridss.pass.decoy.vcf

##########################
## Filter Genotype unknown or ref-ref
##########################



##########################
## Filter Read depth 10 < SV < 1000
##########################



##########################
## Filter size 50 < SV < 10000
##########################


###########################
## Get Coverage and Insert size 
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

# ###########################
# ## Merge SV 
# ###########################

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

############################
### ANNOTATIONS
############################

java -Xmx4g -jar ../bin/snpEff/snpEff.jar -v GRCh37.75 ../data/$1/analysis/all.vcf > ../data/$1/analysis/all.annot.vcf

cat ../data/$1/analysis/all.annot.vcf | grep -v "#" | cut -f8 | cut -d";" -f11 | cut -d"|" -f2 | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//'  | awk -v OFS="\t" '$1=$1' | awk -v FS='\t' -v OFS='\t' '{$2=$2!=""?$2:"NA"}1'  >   ../data/$1/analysis/annotation_counts_no_header 
cat ../data/annotations.header ../data/$1/analysis/annotation_counts_no_header > ../data/$1/analysis/annotation_counts

rm ../data/$1/analysis/annotation_counts_no_header

#########################
## Number of events per sample
#########################

cat ../data/$1/analysis/all.vcf  | grep -v "#"  | cut -f8 | cut -d";" -f4 | cut -d "=" -f2 | sort | uniq -c |   sed -e 's/^[[:space:]]*//'  | awk -v OFS="\t" '$1=$1' | sed "s/^/$1\t/g" > ../data/$1/analysis/number_events_sv_no_header
cat ../data/number_events_sv.header ../data/$1/analysis/number_events_sv_no_header > ../data/$1/analysis/number_events_sv

rm ../data/$1/analysis/number_events_sv_no_header

#########################
## Size distribution of events
#########################

cat ../data/$1/analysis/all.vcf  | grep -v "#"  | cut -f8 | cut -d";" -f3,4 | sed s/SVLEN=//g |  sed s/SVTYPE=//g | sed s/-//g | sed 's/;/\t/g' | sed "s/^/$1\t/g" > ../data/$1/analysis/size_distribution_events_no_header
cat ../data/size_distribution_events.header ../data/$1/analysis/size_distribution_events_no_header > ../data/$1/analysis/size_distribution_events

rm ../data/$1/analysis/size_distribution_events_no_header

#########################
## SV Size Number 
#########################

cat ../data/$1/analysis/size_distribution_events |  tail -n +2 | sort | uniq -c |  sort -nr | sed 's/^ *//g' | sed 's/ /'$'\t''/' > ../data/$1/analysis/sv_size_discovered_no_header
cat ../data/size_count_distribution_events.header ../data/$1/analysis/sv_size_discovered_no_header >  ../data/$1/analysis/sv_size_discovered

rm ../data/$1/analysis/sv_size_discovered_no_header

#########################
## Number of events per tool
#########################

cat ../data/$1/delly_out/delly.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_delly
cat ../data/$1/lumpy_out/lumpy.vcf  | grep -v "#" |  wc -l > ../data/$1/analysis/number_events_lumpy
cat ../data/$1/manta_out/manta.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_manta
cat ../data/$1/gridss_out/gridss.vcf  | grep -v "#" | wc -l > ../data/$1/analysis/number_events_gridss

#########################
## Event Intersection
#########################

#cat ../data/$1/analysis/4.vcf | TODO


# #delly
# cat ../data/$1/delly_out/delly.vcf  | grep -v "#" | cut -f1,2 > ../data/delly_chr_start
# cat ../data/$1/delly_out/delly.vcf | grep -v "#" | cut -f8 | cut -d";" -f5 | cut -d"=" -f2 > ../data/delly_end
# paste ../data/delly_chr_start ../data/delly_end > ../data/delly.bed


# #manta
# cat ../data/$1/manta_out/manta.vcf  | grep -v "#" |  cut -f1,2> ../data/manta_chr_start
# cat ../data/$1/manta_out/manta.vcf | grep -v "#" | cut -f8 | cut -d";" -f1 | cut -d"=" -f2 > ../data/manta_end
# paste ../data/manta_chr_start ../data/manta_end > ../data/manta.bed


# #lumpy 
# cat ../data/$1/lumpy_out/lumpy.vcf  | grep -v "#" | grep -v "BND" |  cut -f1,2 > ../data/lumpy_chr_start
# cat ../data/$1/lumpy_out/lumpy.vcf  | grep -v "#" | grep -v "BND" | cut -f8 | cut -d";" -f4 | cut -d"=" -f2 > ../data/lumpy_end
# paste ../data/lumpy_chr_start ../data/lumpy_end > ../data/lumpy.bed


# #gridss
# cat ../data/$1/gridss_out/gridss.vcf  | grep -v "#" |  cut -f1,2 > ../data/gridss_chr_start
# cat ../data/$1/gridss_out/gridss.vcf | grep -v "#" | > ../data/gridss_end
# paste ../data/gridss_chr_start ../data/gridss_end > ../data/gridss.bed

# # calculate intersections

# bedtools intersect -wa -wb -a delly.bed -b manta.bed lumpy.bed -filenames |  cut -f4 | cut -d"." -f1 | sort | uniq -c
# bedtools intersect -wa -wb -a manta.bed -b delly.bed lumpy.bed -filenames |  cut -f4 | cut -d"." -f1 | sort | uniq -c
# bedtools intersect -wa -wb -a lumpy.bed -b delly.bed manta.bed -filenames |  cut -f4 | cut -d"." -f1 | sort | uniq -c


#########################
## Strucutral Variants per Read-Depth
#########################

# #sv
# cat ../data/$1/delly_out/delly.vcf | grep -v "#"| cut -f8 | cut -d";" -f2  | cut -d"=" -f2 > ../data/$1/delly_out/delly.vcf/analysis/delly_sv
# #size # TODO
# cat ../data/$1/delly_out/delly.vcf | grep -v "#"| cut -f8 | cut -d";" -f2  | cut -d"=" -f2 > ../data/$1/delly_out/delly.vcf/analysis/delly_sv
# #PE
# cat ../data/$1/delly_out/delly.vcf | grep -v "#"| cut -f8 | cut -d";" -f6  | cut -d"=" -f2 > ../data/$1/delly_out/delly.vcf/analysis/delly_PE
# #SR
# cat ../data/$1/delly_out/delly.vcf | grep -v "#"| cut -f8 | cut -d";" -f13 | cut -d"=" -f2 > ../data/$1/delly_out/delly.vcf/analysis/delly_SR

# #calculate sum SV RD
# paste delly_PE delly_SR |  awk -v FS='\t' -v OFS='\t' '{$2=$2!=""?$2:"0"}1' | awk 'NR==1 { print "sum"; next } { print $1 + $2 }' > ./data/$1/delly_out/delly.vcf/analysis/sum_sv_rd

# paste ../data/$1/delly_out/delly.vcf/analysis/delly_sv ../data/$1/delly_out/delly.vcf/analysis/sum_sv_rd | sed 's/^/NA12878\t/g' > ../data/$1/delly_out/delly.vcf/analysis/sv_rd

# # count events 

# # Add Header
# sed -i '1iSAMPLE SV SIZE RD' ../data/$1/analysis/sv_rd


# ############################
# ## CATALOG OF VARIANTS (RARE)
# ############################

## create ALL BED ##

# ../bin/SURVIVOR/Debug/./SURVIVOR vcftobed ../data/$1/analysis/all.annot.vcf  0 -1 ../data/$1/all.annot.bed

## GET SV-GENOMAD VCF ##

#cat ../data/$1/analysis/all.annot.vcf | grep -v "#" | cut -f1,2 > ../data/$1/analysis/all_chr_start
#cat ../data/$1/analysis/all.annot.vcf | grep -v "#" | cut -f8 | cut -d";" -f7 | cut -d"=" -f2 > ../data/$1/analysis/all_end
#paste ../data/$1/analysis/all_chr_start ../data/$1/analysis/all_end > all.bed

## PULL GNOMAD DB  https://gnomad.broadinstitute.org/downloads

#bedtools intersect -wa -wb -a ../data/$1/all.annot.bed -b ../data/gnomad_v2.1_sv.sites.bed -filenames

#########################
## GENERATE R MARKDOWN
#########################

sed  "s/sample_template/$1/g" < sample_markdown_template.rmd > $1.rmd

####################
## GENERATE HTML FROM MARKDOWN
#####################

Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); rmarkdown::render('"$1".rmd')"



