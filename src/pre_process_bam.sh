#!/bin/bash

##############################################################
##############################################################
# Configuration

BQSR="FALSE"

## Not applied
#POST_BQSR_ANALYSIS="FALSE"
#HQ_FILTER_STEP="FALSE"
#HQ_THRESHOLD=20

##############################################################

filename="$(basename -- "${1%.*}")"
echo "The sample name is '$filename'."

TIME=$(date +%Y%m%d%H%M%S)
enddir=${filename}_${TIME}
mkdir -p $enddir

path=${1}
fullname=$(echo $(readlink -f "$path"))

##############################################################
##############################################################

### INFORMATION OF THE BAM ###

####################
## Get Mapper
####################

echo "Start get mapper"

# Mapper tool
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep '^@PG' | grep -wo bwa | head -n1 > ${enddir}/${filename}.pg_tags

# Mapper method
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep '^@PG' | grep bwa | grep -woE 'aln|mem' | head -n1 >> ${enddir}/${filename}.pg_tags

# Mapper version
/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep '^@PG' | grep bwa | grep -E 'aln|mem' | grep -o "VN:.*" | awk '{print $1}' | cut -c4- | head -n1 >> ${enddir}/${filename}.pg_tags

echo "End get mapper"

###################
### Get read groups
###################

echo "Start get RG"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep '^@RG' | grep -o "ID:.*" | cut -f1 > ${enddir}/${filename}.rg_tags

echo "End get RG"

################
### Get quality threshold
################

echo "Start get quality threshold"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep '^@PG' | grep bwa | grep -o "\-q.*" | cut -c4- | awk '{print $1}' | head -n1 > ${enddir}/${filename}.quality_threshold

echo "End get quality threshold"

###############################################################
###############################################################

# START OF THE PIPELINE


# Getting if the bam file is sorted by coordinates and if its already marked by duplicates

SORT_COORDINATE=$(/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname} | grep "^@HD" | grep "coordinate")

MARK_DUPLICATES=$(/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -f 0x400 ${fullname} | head -n1) 

if [ ! -z "$SORT_COORDINATE" ];
then

#####################
### Sort by coordinate (samtools)
#####################

    echo "Start sort by coordinate (samtools)"

    /iso/tmp/samtoolsINST-1.3.1/bin/samtools sort -m 1500M --threads 35 -T /vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/temp_files ${fullname} > ${enddir}/${filename}.sort.bam

    echo "End sort by coordinate (samtools)"

else
    echo "¡ASSUMING THAT IS ORDERED BY COORDINATE!"
    echo "If not, you need to set the sort configuration parameter to 'TRUE'"
    cp ${fullname} ${enddir}/${filename}.sort.bam
fi

###############################################################

### GATK BEST PRACTICES ###

if [ ! -z "$MARK_DUPLICATES" ];
then

#####################
### Mark Duplicates (bam file ordered by coordinate) and remove duplicates. Use MarkDuplicatesSpark if you want to run it in parallel.
#####################

    echo "Start mark duplicates"

    /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk MarkDuplicates -I ${enddir}/${filename}.sort.bam -O ${enddir}/${filename}.sort.marked_duplicates.bam -M ${enddir}/${filename}.sort.marked_duplicates_metrics.txt --TMP_DIR /vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/temp_files --VALIDATION_STRINGENCY LENIENT

    # Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

    echo "End mark duplicates"

else
    echo "¡ASSUMING THAT THE FILE IS MARKED BY DUPLICATES!"
    echo "If not, you need to set the sort configuration parameter to 'TRUE'"
    cp ${enddir}/${filename}.sort.bam ${enddir}/${filename}.sort.marked_duplicates.bam

fi
    
################################################################

if [[ "$BQSR" -eq "TRUE" ]] ;
then

######################
### BaseRecalibrator before
######################

    echo "Start BaseRecalibrator before"

    /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk BaseRecalibrator -I ${enddir}/${filename}.sort.marked_duplicates.bam -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta --known-sites /vault/bio_scratch/arnau/structural_variants/gnomad_v2.1_sv.sites.vcf.gz -O ${enddir}/${filename}.recal_data.table1

    echo "End BaseREcalibrator before"

#####################
### ApplyBQSR
#####################

    echo "Start ApplyBQSR"

    /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk ApplyBQSR -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta -I ${enddir}/${filename}.sort.marked_duplicates.bam --bqsr-recal-file ${enddir}/${filename}.recal_data.table1 -O ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam  --add-output-sam-program-record true 

#--create-output-bam-index true #--preserve-qscores-less-than 10 
    echo "End ApplyBQSR"

else
	echo "¡BQSR IS NOT APPLIED!"
fi

###############################################################

#####################
## Remove intermediate files
#####################

rm ${enddir}/${filename}.sort.bam

if [ "$BQSR" = "TRUE" ];
then
    rm ${enddir}/${filename}.sort.marked_duplicates.bam
fi

#rm ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam

################################################################

#######################
## Create Index
#######################

echo "Start create index"

/iso/tmp/samtoolsINST-1.3.1/bin/samtools index ${enddir}/${filename}.sort.marked_duplicates*bam

echo "End create index"

###############################################################

# FINAL STATISTICS

cd ${enddir}

if [ "$BQSR" = "TRUE" ];
then

    bash ../calculate_statistics.sh ${enddir}/${filename}.sort.marked_duplicates.BQSR

    bash ../calculate_depth_coverage.sh ${enddir}/${filename}.sort.marked_duplicates.BQSR

else

    bash ../calculate_statistics.sh ${enddir}/${filename}.sort.marked_duplicates

    bash ../calculate_depth_coverage.sh ${enddir}/${filename}.sort.marked_duplicates
fi

cd ..

#####################
# How many records are in the filtered BAM compared to the original? How many read pairs does this represent?
#####################

echo "How many records are in the filtered BAM compared to the original:" > ${enddir}/${filename}.records_last
echo "`/iso/tmp/samtoolsINST-1.3.1/bin/samtools view -c ${fullname}.bam` `/iso/tmp/samtoolsINST-1.3.1/bin/samtools view -c ${enddir}/${filename}.sort.marked_duplicates*.bam`" | awk '{printf "%0.2f%%\n", 100*$2/$1}' >> ${enddir}/${filename}.records_last

###############################################################
###############################################################

### Optionals steps

####################
# Filter by Quality > 20 (-q)
####################

#echo 'Sart filter reads by quality' && /iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -b -q 20 ${enddir}/${filename}.sort.marked_duplicates*bam > ${enddir}/${filename}.sort.marked_duplicates*HQ.bam && echo 'End filter reads'

####################################

### Analysis of the recalibration step

####################
### BaseRecalibrator after
######################

#echo "Start BaseRecalibrator after" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk BaseRecalibrator -I ${enddir}/${filename}.sort.marked_duplicates*HQ.bam -R /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta --known-sites /vault/bio_scratch/arnau/structural_variants/gnomad_v2.1_sv.sites.vcf.gz -O ${enddir}/${filename}.recal_data.table2 && echo "End BaseREcalibrator after"

######################
### AnalyzeCovariates
#####################

#echo "Start AnalyzeCovariates" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk AnalyzeCovariates -before ${enddir}/${filename}.recal_data.table1 -after ${enddir}/${filename}.recal_data.table2 -plots ${enddir}/${filename}.AnalyzeCovariates_last.pdf && echo "End AnalyzeCovariates"

###############################################################
###############################################################

rm /vault/mauricio/bio_team/SV/sv-callers/workflow/data/bam/temp_files/* 

echo "Check the results in '${enddir}' folder!!"
