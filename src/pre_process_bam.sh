#!/bin/bash

# INFORMATION OF THE BAM

####################
## Get Mapper
####################

echo "Start get mapper"

# Mapper tool
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H $1.bam | grep '^@PG' | grep -wo bwa | head -n1 > $1.pg_tags

# Mapper method
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H $1.bam | grep '^@PG' | grep bwa | grep -woE 'aln|mem' | head -n1 >> $1.pg_tags

# Mapper version
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H $1.bam | grep '^@PG' | grep bwa | grep -E 'aln|mem' | grep -o "VN:.*" | awk '{print $1}' | cut -c4- | head -n1 >> $1.pg_tags

echo "End get mapper"

###################
### Get read groups
###################

echo "Start get RG"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H $1.bam | grep '^@RG' | grep -o "ID:.*" | cut -f1 > $1.rg_tags

echo "End get RG"

################
### Get quality threshold
################

echo "Start get quality threshold"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H $1.bam | grep '^@PG' | grep bwa | grep -o "\-q.*" | cut -c4- | awk '{print $1}' | head -n1 > $1.quality_threshold

echo "End get quality threshold"

###############################################################
###############################################################

# START OF THE PIPELINE

#####################
### Sort by coordinate (samtools)
#####################

echo "Start sort by coordinate (samtools)"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools sort -m 1500M --threads 35 -T /vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/temp_files $1.bam > $1.sort.bam

echo "End sort by coordinate (samtools)"

###############################################################

### GATK BEST PRACTICES ###

#####################
### Mark Duplicates (bam file ordered by coordinate) and remove duplicates. Use MarkDuplicatesSpark if you want to run it in parallel.
#####################

echo "Start mark duplicates"

/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk MarkDuplicates -I $1.sort.bam -O $1.sort.marked_duplicates.bam -M $1.sort.marked_duplicates_metrics.txt --TMP_DIR /vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/temp_files --VALIDATION_STRINGENCY LENIENT

# Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

echo "End mark duplicates"

######################
### BaseRecalibrator before
######################

echo "Start BaseRecalibrator before"

/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk BaseRecalibrator -I $1.sort.marked_duplicates.bam -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta --known-sites /vault/arnau/structural_variants/gnomad_v2.1_sv.sites.vcf.gz -O $1.recal_data.table1

echo "End BaseREcalibrator before"

#####################
### ApplyBQSR
#####################

echo "Start ApplyBQSR"

/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk ApplyBQSR -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta -I $1.sort.marked_duplicates.bam --bqsr-recal-file $1.recal_data.table1 -O $1.sort.marked_duplicates.BQSR.bam  --add-output-sam-program-record true #--create-output-bam-index true #--preserve-qscores-less-than 10 

echo "End ApplyBQSR"

#######################
## Create Index
#######################

echo "Start create index"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools index $1.sort.marked_duplicates.BQSR.bam

echo "End create index"

###############################################################

# FINAL STATISTICS

######################
## Calculate Statistics
######################

bash calculate_statistics.sh $1.sort.marked_duplicates.BQSR

######################
## Calculate Depth of Coverage
######################

bash calculate_depth_coverage.sh $1.sort.marked_duplicates.BQSR

#####################
# How many records are in the filtered BAM compared to the original? How many read pairs does this represent?
#####################

echo "How many records are in the filtered BAM compared to the original:" > $1.records_last
echo "`/iso/tmp/samtoolsINST-1.3.1/bin/samtools view -c $1.bam` `/iso/tmp/samtoolsINST-1.3.1/bin/samtools view -c $1.sort.marked_duplicates.BQSR.bam`" | awk '{printf "%0.2f%%\n", 100*$2/$1}' >> $1.records_last

###############################################################
###############################################################

### Optionals steps

####################
# Filter by Quality > 20 (-q)
####################

#echo 'Sart filter reads by quality' && /iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -b -q 20 $1.sort.marked_duplicates.BQSR.bam > $1.sort.marked_duplicates.BQSR.HQ.bam && echo 'End filter reads'

####################################

### Analysis of the recalibration step

####################
### BaseRecalibrator after
######################

#echo "Start BaseRecalibrator after" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk BaseRecalibrator -I $1.sort.marked_duplicates.BQSR.HQ.bam -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta --known-sites /vault/arnau/structural_variants/gnomad_v2.1_sv.sites.vcf.gz -O $1.recal_data.table2 && echo "End BaseREcalibrator after"

######################
### AnalyzeCovariates
#####################

#echo "Start AnalyzeCovariates" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk AnalyzeCovariates -before $1.recal_data.table1 -after $1.recal_data.table2 -plots $1.AnalyzeCovariates_last.pdf && echo "End AnalyzeCovariates"

###############################################################
###############################################################

#####################
## Remove intermediate files
#####################

rm $1.sort.bam
rm $1.sort.marked_duplicates.bam
#rm $1.sort.marked_duplicates.BQSR.bam


