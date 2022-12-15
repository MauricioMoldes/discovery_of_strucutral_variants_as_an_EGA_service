#!/usr/bin/env bash

#####################
## Pull data
#####################

bash pull_data.sh $1

######################
## Post-Processing data
#######################

bash post_processing.sh $1

###########################
## Coverage and Insert size 
###########################

bash coverage_insert_size.sh $1

###########################
## Merge SV 
###########################

bash merge_sv.sh $1

############################
### Annotations
############################

bash annotations.sh $1

#########################
## Number of events per sample
#########################

bash number_events.sh $1

#########################
## Size distribution of events
#########################

bash size_distribution_events.sh $1

#########################
## SV Size Number 
#########################

bash sv_size_number.sh $1

#########################
## Number of events per tool
#########################

bash number_events_tool.sh $1

#########################
## GENERATE R MARKDOWN
#########################

sed  "s/sample_template/$1/g" < sample_markdown_template.rmd > ../results/$1.rmd

####################
## GENERATE HTML FROM MARKDOWN
#####################

cd ../results/

Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); rmarkdown::render('"$1".rmd')"



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



