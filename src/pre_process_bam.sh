#!/bin/bash

##############################################################

## STRUCTURAL VARIANTS PRE-PROCESS PIPELINE ##

Author="Arnau Soler Costa"
Creation_date="09/02/2023"
Update_date="02/03/2023"
Version="v0.0.1"

##############################################################

# Script to pre-process the bam files. 
# This pipeline follows the GATK (gatk-4.2.6.0) best practices. The execution of this script allows to obtain a final bam file: sorted (samtools v.1.3.1), marked the duplicates (GATK) and perform an optional step called 'Base Quality Score Recalibration' (BQSR from GATK). For the 'BQSR' step we are using the information on gnomAD (Broad Institute) structural variants vcf as the known sites of real variants (https://gnomad.broadinstitute.org/downloads/#v2-structural-variants).

# The pipeline has 'logic coding' skipping the steps of sorting and marking the duplicates if the initial bam file is already sorted and marked. Moreover, it's possible to select if you want to sort and/or mark the duplicates again, or force the script to do it if there is a problem with the initial bam. It's possible to choose if the 'BQSR' step is run or not, and when running this step, it's possible to choose if removing all intermidiate files or keeping the bam file sorted and marked the duplicates while the one sorted, marked the duplicates and 'BQSR', among other parameters.

# The script has a set of commands to return a set of files with information and result statistics of the initial and final bams. It's also possible to choose to no return these final informative files. 

# The script creates a new and final directory where is possible to find the all the results.

#  No read is trimmed or removed in this pipeline.

##############################################################
##############################################################

###################
# Configuration variables and errors
###################

# Default variables values

BAM_FILE="" # Path to input bam file
BQSR=false # True = Execute BQSR function, False = no execution BQSR.
FORCE_SORT=true # True = Force sort bam file by coordinates.
FORCE_MARK_DUPLICATES=true # True = Force mark duplicates the bam file.
RM_SORT_MARKED_FILE=false # True = Remove the intermidate file 'sort.marked_duplicates'. Only set 'True' if the 'BQSR' is also set True.
INITIAL_INFO=true # True = Get initial information of the initial bam file (True).
FINAL_STATISTICS=true # Get final statistics of the final bam file, and comparison with the initial bam file (True).
TMP_DIR="" # Set the path to the temprorary files directory
RM_TMP_FILES=true # Set 'true' to remove the temporary files. False no remove.
INDEX_INPUT=false # Set 'true' to index the input bam file 
REFERENCE_GENOME="" # Path to the reference genome of the input bam file
CLEAN_SAM=false # To execute the CleanSam fuction (true). Default false.
VALIDATE_SAM=false # To execute the ValidateSamFile fuction (true). Default false.

## Not applied
#POST_BQSR_ANALYSIS="FALSE"
#HQ_FILTER_STEP="FALSE"
#HQ_THRESHOLD=20

##############################################################

# Control errors (Add -e to exit with non-zero errors. Add -x to debug)
set -uoh pipefail

# Sets the locale environment variable to "C" in the Bash shell
export LC_ALL=C

##############################################################

# Function to print a help message.
usage() {                                 
  echo "Usage: $0 
            
            Required parameters:

            -I/--INPUT_BAM_FILE --- Path to the bam to be processed (ie. filename path).
	    
	    -R/--REFERENCE_GENOME --- Path to reference genome of the input bam file (ie. b37, hg19, GRCh38,...).
            
            -T/--TMP_DIR --- Path to the directory for the temporary files. (ie. directory path). No default. 
 
	    Optional parameters:

	    -Q/--BQSR --- If TRUE, BQSR will be performed (ie. true/false). Default 'false'.
	    
	    -S/--FORCE_SORT_COORDINATES --- Force sort the bam file by coordinates (ie. true/false). Default 'true'.
	    
	    -M/--FORCE_MARK_DUPLICATES --- Force mark the duplicates of the bam file. Default 'true'. 
	    
	    -r/--rm_sort_marked_file --- Remove the sorted and marked bam file when you perform the BQSR as a final step (ie. true/false). Default 'false'.
	    
	    -i/--initial_info --- Show the initial information of the input bam file. Default true. (ie. true/false). Default 'true'.
	    
	    -f/--final_statistics --- Perform and show some final statistics from the last obtained bam file and comparisons with the input bam file. Default true. (ie. true/false). Default 'true'.
	    
	    -t/--rm_tmp_files --- Delete the temporary files. Default true. (ie. true/false). Default 'true'.
	    
	    -x/--index_input_file --- If TRUE, indexing of the input bam file will be performed (ie. true/false). Default 'false'.

	    -C/--CLEAN_SAM --- To clean the SAM/BAM file by the function CleanSam (Picard - GATK) (ie. true/false). Default 'false'.
	    
	    -V/--VALIDATE_SAM_FILE --- To validate the SAM/BAM file by the function ValidateSamFile (Picard - GATK) (ie. true/false). Default 'false'.

	    -h/--help --- Show usage of the script
	    
	    -v/--version --- Show version of the script" 1>&2
}

# Function: Exit with error
exit_abnormal() {           
  usage
  exit 1
}

# Control variables
while [[ "$#" -gt 0 ]]; do    
                              
  case $1 in                    
    -I|--INPUT_BAM_FILE) BAM_FILE="$2"; 
      
      # Check if the path exists
      if ! test -f "$BAM_FILE"; then
        
	  # Check if the file is a '.bam'
          if [ "${fullname: -4}" = ".bam" ];
          then

              header=$(/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 -H ${fullname})
              #body=$(/iso/tmp/samtoolsINST-1.3.1/bin/samtools view --threads 35 ${fullname} | head -n1)

              if [ -z "${header}" ];
              then
                  echo "ERROR: File does not have a valid header and/or body"
                  exit_abnormal
              else
                  echo "Valid 'bam' file"
              fi
          else
              echo "ERROR: File does not have the proper extension (.bam) or it can't be read"
              exit_abnormal
          fi

	     echo "--- CAN'T START THE PIPELINE, FILE OR PATH ARE NOT CORRECT ---"
             exit_abnormal
      fi
      shift
      ;;  
    -Q|--BQSR) BQSR="$2";                     
      if ! [[ "${BQSR,,}" == "true" ]] && ! [[ "${BQSR,,}" == "false" ]]; then   # if $BQSR is not boolean:
	echo "ERROR: BQSR (-Q) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -S|--FORCE_SORT_COORDINATES) FORCE_SORT="$2";                     
      if ! [[ "${FORCE_SORT,,}" == "true" ]] && ! [[ "${FORCE_SORT,,}" == "false" ]]; then   # if $FORCE_SORT is not boolean:
	echo "ERROR: FORCE_SORT_COORDINATE (-S) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -M|--FORCE_MARK_DUPLICATES) FORCE_MARK_DUPLICATES="$2";
      if ! [[ "${FORCE_MARK_DUPLICATES,,}" == "true" ]] && ! [[ "${FORCE_MARK_DUPLICATES,,}" == "false" ]]; then   # if $FORCE_MARK_DUPLICATES is not boolean:
        echo "ERROR: FORCE_MARK_DUPLICATES (-M) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -r|--rm_sort_marked_file) RM_SORT_MARKED_FILE="$2";                     # Set $RM_SORT_MARKED_FILE to specified value.
      if ! [[ "${RM_SORT_MARKED_FILE,,}" == "true" ]] && ! [[ "${RM_SORT_MARKED_FILE,,}" == "false" ]]; then   # if $RM_SORT_MARKED_FILE is not boolean:
        echo "ERROR: rm_sort_marked_file (-r) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -i|--initial_info) INITIAL_INFO="$2";                     # Set $INITIAL_INFO to specified value.
      if ! [[ "${INITIAL_INFO,,}" == "true" ]] && ! [[ "${INITIAL_INFO,,}" == "false" ]]; then   # if $INITIAL_INFO is not boolean:
        echo "ERROR: initial_info (-i) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -f|--final_statistics) FINAL_STATISTICS="$2";                     # Set $FINAL_STATISTICS to specified value.
      if ! [[ "${FINAL_STATISTICS,,}" == "true" ]] && ! [[ "${FINAL_STATISTICS,,}" == "false" ]]; then   # if $FINAL_STATISTICS is not boolean:
        echo "ERROR: final_statistics (-f) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -T|--TMP_DIR) TMP_DIR="$2";
      if ! [ -d "$TMP_DIR" ]; then
        echo "ERROR: Wrong path to directory of temporary files"
        exit_abnormal
      fi
      shift
      ;;   
    -t|--rm_tmp_files) RM_TMP_FILES="$2";                     # Set $RM_TMP_FILES to specified value.
      if ! [[ "${FORCE_SORT,,}" == "true" ]] && ! [[ "${FORCE_SORT,,}" == "false" ]]; then   # if $RM_TMP_FILES is not boolean:
        echo "ERROR: RM_TMP_FILES (-t) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -x|--index_input_file) INDEX_INPUT="$2";
      if ! [[ "${INDEX_INPUT,,}" == "true" ]] && ! [[ "${INDEX_INPUT,,}" == "false" ]]; then   
	      echo "ERROR: index_input_file(-x) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -R|--REFERENCE_GENOME) REFERENCE_GENOME="$2";
          
      # Check if the path exists
      if ! test -f "$REFERENCE_GENOME"; then
	  echo "ERROR: --REFERENCE_GENOME (-R)  path is not correct."
	  exit_abnormal
      fi
      shift
      ;;
    -C|--CLEAN_SAM) CLEAN_SAM="$2";
      if ! [[ "${CLEAN_SAM,,}" == "true" ]] && ! [[ "${CLEAN_SAM,,}" == "false" ]]; then
	echo "ERROR: CLEAN_SAM (-C) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -V|--VALIDATE_SAM_FILE) VALIDATE_SAM="$2";
      if ! [[ "${VALIDATE_SAM,,}" == "true" ]] && ! [[ "${VALIDATE_SAM,,}" == "false" ]]; then
        echo "ERROR: VALIDATE_SAM_FILE (-V) must be a boolean, true/false."
        exit_abnormal
      fi
      shift
      ;;
    -h|--help)
      exit_abnormal
      ;;
    -v|--version)
      echo "Author: ${Author}"
      echo "Creation date: ${Creation_date}"
      echo "Update date: ${Update_date}"
      echo "Version: ${Version}"
      exit 0
      ;;    
    :)
      echo "ERROR: The parameter $1 is wrong";
      exit_abnormal
      ;;
    *) # If unknown (any other) option:
      echo "ERROR: One or more parameters are wrong.";
      exit_abnormal # Exit abnormally.
      ;;
  esac
  shift
  
done

# Check if input, reference genomes and temporary files directory are specified (those are mandatory parameters).

# Check if a specific parameter was passed to the script
if [ -z "$BAM_FILE" ]; then
    echo "ERROR: Parameter -I/--INPUT_BAM_FILE is required."
    exit_abnormal
fi

# Check if a specific parameter was passed to the script
if [ -z "$REFERENCE_GENOME" ]; then
    echo "ERROR: Parameter -R/--REFERENCE_GENOME is required."
    exit_abnormal
fi

# Check if a specific parameter was passed to the script
if [ -z "$TMP_DIR" ]; then
    echo "ERROR: Parameter -T/--TMP_DIR is required."
    exit_abnormal
fi

##################################################################

# Set names and path directories

# Need to be 'sudo' to run the script
#[[ $(id -u) -eq 0 ]] || { echo >&2 "Must be root to run the script properly"; exit 1; }

# Names, filenames and paths
# Basename of the file
filename=$(basename $BAM_FILE .bam)

# Full path to file
fullname=$(echo $(readlink -f "$BAM_FILE"))

# Create results folder (end directory)
TIME=$(date +%Y%m%d%H%M%S)
enddir=${filename}_${TIME}
mkdir -p $enddir

# Create error logs file
touch ${enddir}/${filename}_split_sample.log

###############################################################

echo "##### START PIPELINE #####"

####################
## Check multi-sample bam file
####################

# Finf if the bam is multi sample

multi_sample=$(samtools view -H ${fullname} | grep "^@RG" | grep "SM" | sed 's/^.*\(SM:.*\).*$/\1/' | cut -f1 | sort | uniq | wc -l)

# If it's multi-sample bam, split it
if [[ $multi_sample -gt 1 ]];then
    
    echo "# SPLIT SAMPLE LOGS FILE FOR ${filename}.bam" &>> ${enddir}/${filename}_split_sample.log	
    echo "The input file is a multi-sample bam file" &>> ${enddir}/${filename}_split_sample.log
    echo "Splitting the bam file by samples" &>> ${enddir}/${filename}_split_sample.log
    
    # Split reads function
    gatk SplitReads -I ${fullname} -O ${enddir} --split-sample true --create-output-bam-index false --read-validation-stringency LENIENT --tmp-dir ${TMP_DIR} --add-output-sam-program-record true &>> ${enddir}/${filename}_split_sample.log
    
    # List of samples
    find ${enddir} -type f -name "*bam" > ${enddir}/list_of_sample.txt

else
    #echo "The input file is a single sample bam file"
    echo ${fullname} > ${enddir}/list_of_sample.txt	
fi

# Remove split_sample.log file if it's empty (single sample input BAM file)
if ! [ -s "${enddir}/${filename}_split_sample.log" ]; then
  rm ${enddir}/${filename}_split_sample.log
fi

###############################################################
###############################################################

## Start of pipeline

for sample in `cat ${enddir}/list_of_sample.txt`
do

# Names, filenames and paths

# Basename of the file
filename=$(basename $sample .bam)

# Full path to file
fullname=$(echo $(readlink -f "$sample"))

# Index file name
INDEX="${fullname}.bai"

# Create error logs file
touch ${enddir}/${filename}.log

# Config parameters for logs (info. https://serverfault.com/questions/103501/how-can-i-fully-log-all-bash-scripts-actions)
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>${enddir}/${filename}.log 2>&1
echo -e "\n LOGS FILE FOR ${filename}.bam \n"
# Everything below will go to the file '$sample.log':
# Only work inside the loop (because is inside it)
# Adding '>&3' to the end of the command, the output goes to the console, not the 'log' file 
# Other way: Adding '&>> ${enddir}/${filename}.log' to the end of each command

########################

echo "## PROCESSING FILE: ${filename}.bam" >&3 # '>&3' to show to the console

echo "# PROCESSING FILE: ${filename}.bam" # Print to the log file

# Clean sam (test)
if [ "${CLEAN_SAM,,}" = "true" ];
then
    echo "Start CleanSam"
     
    gatk CleanSam --INPUT ${fullname} --OUTPUT ${enddir}/${filename}.bam --VALIDATION_STRINGENCY LENIENT --TMP_DIR ${TMP_DIR}
    
    echo "CleanSam done"

fi

# Validate Sam file before cleaning
if [ "${VALIDATE_SAM,,}" = "true" ];
then
    #echo "Start ValidateSamFile - MODE=SUMMARY"
    #/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk ValidateSamFile -I ${fullname} --MODE SUMMARY --TMP_DIR ${TMP_DIR
    
    echo "Example provided to perform the Validation of the BAM file by the user using ValidateSamFile (Picard - GATK) tool and be able to fix the errors/warnings if needed:
          
          SUMMARY MODE:

                  ValidateSamFile -I input.bam --MODE SUMMARY

          VERBOSE MODE (seeing only errors):

	          ValidateSamFile -I input.bam --MODE VERBOSE --IGNORE_WARNINGS true
       
	  VERBOSE MODE (seeing only warnings):

                  ValidateSamFile -I input.bam --MODE VERBOSE --IGNORE type
	
	  https://gatk.broadinstitute.org/hc/en-us/articles/360036854731-ValidateSamFile-Picard-
          https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/solutions-to-problems/7571-Errors-in-SAMBAM-files-can-be-diagnosed-with-ValidateSamFile" 

fi

# Index if necessary
# Index input bam file
# Index file if its not already indexed
if test -f "${INDEX}"; then
    echo "${filename}.bam.bai already exists"
else
    if [ "${INDEX_INPUT,,}" = "true" ];
    then

        echo "Start indexing ${filename}.bam" 

        samtools index ${fullname} 

        echo "End index"
    else
	echo "WARNING: No indexed file of input bam in the folder"
    fi

fi

##############################################################
##############################################################

### INFORMATION OF THE BAM ###

if [ "${INITIAL_INFO,,}" = "true" ];
then
    
    echo "## START INITIAL INFORMATION" 

    # Get all PG tags    
    samtools view -H ${fullname} | grep '^@PG' | head -n1 > ${enddir}/${filename}.pg_tags
     
    # Get Mapper
    # Mapper tool
    echo "BWA mapper" > ${enddir}/${filename}.mapper
    
    # Mapper
    samtools view -H ${fullname} | grep '^@PG' | grep -wo bwa | head -n1 >> ${enddir}/${filename}.mapper

    # Mapper method
    samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -woE 'aln|mem' | head -n1 >> ${enddir}/${filename}.mapper

    # Mapper version
    samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -E 'aln|mem' | grep -o "VN:.*" | awk '{print $1}' | cut -c4- | head -n1 >> ${enddir}/${filename}.mapper

    # Get read groups
    samtools view -H ${fullname} | grep '^@RG' | grep -o "ID:.*" | cut -f1 > ${enddir}/${filename}.rg_tags

    # Get quality threshold
    samtools view -H ${fullname} | grep '^@PG' | grep bwa | grep -o "\-q.*" | cut -c4- | awk '{print $1}' | head -n1 > ${enddir}/${filename}.quality_threshold

    echo "## INITIAL INFORMATION DONE" 

fi

###############################################################
###############################################################

# START OF THE PIPELINE

echo "### START OF THE PIPELINE" 

#####################
### Sort by coordinate (samtools)
#####################

# Getting if the bam file is sorted by coordinates

SORT_COORDINATE=$(timeout 5 samtools view -H ${fullname} | grep "^@HD" | grep "coordinate")

if [ -z "$SORT_COORDINATE" ] || [ "${FORCE_SORT,,}" = "true" ];
then

    echo "## START SORT BY COORDINATES" 
    
    # SAMTOOLS SORT
    #samtools sort -m 1500M --threads 35 -T ${TMP_DIR} ${fullname} --output-fmt BAM -o ${enddir}/${filename}.sort.bam
    samtools sort -T ${TMP_DIR} ${fullname} --output-fmt BAM -o ${enddir}/${filename}.sort.bam

    
    # GATK SORT
    #/vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk SortSam INPUT=${fullname} OUTPUT=${enddir}/${filename}sort.bam SORT_ORDER=coordinate --TMP_DIR /vault/bio-scratch/arnau/structural_variants/tmp_files --VALIDATION_STRINGENCY LENIENT 
    
    echo "## SORT BY COORDINATES DONE" 

else
    echo "WARNING: File is already sorted by coordinates. Samtools sort is not applied" 
    
    echo "Copying file to the correct folder" 
    cp ${fullname} ${enddir}/${filename}.sort.bam
    echo "File is copyied correctly" 

fi

###############################################################

#####################
### Mark Duplicates (bam file ordered by coordinate) and remove duplicates. Use MarkDuplicatesSpark if you want to run it in parallel.
#####################

# Getting if the bam file is marked by duplicates 

MARK_DUPLICATES=$(timeout 15 samtools view -f 0x400 ${enddir}/${filename}.sort.bam | head -n1)

if [ -z "$MARK_DUPLICATES" ] || [ "${FORCE_MARK_DUPLICATES,,}" = "true" ];
then

#####################
### Mark Duplicates (bam file ordered by coordinate) and remove duplicates. Use MarkDuplicatesSpark if you want to run it in parallel.
#####################

    echo "## START MARK DUPLICATES" 

    gatk MarkDuplicates -I ${enddir}/${filename}.sort.bam -O ${enddir}/${filename}.sort.marked_duplicates.bam -M ${enddir}/${filename}.sort.marked_duplicates_metrics.txt --TMP_DIR ${TMP_DIR} --VALIDATION_STRINGENCY LENIENT --ADD_PG_TAG_TO_READS true --ASSUME_SORTED true 

    # Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. For memory issues (ie. --java-options "-Xmx5g -Djava.io.tmpdir=${TMP_DIR}")

    echo "## MARK DUPLICATES DONE" 

else

    echo "WARNING: File is already marked by duplicates. MarkDuplicates (GATK) is not applied" 
    
    echo "Copying file to the correct folder" 
    cp ${enddir}/${filename}.sort.bam ${enddir}/${filename}.sort.marked_duplicates.bam
    echo "File is copyied correctly" 

fi
    
################################################################

######################
### BaseRecalibrator before
######################

## /vault/mauricio/bio_team/SV/sv-callers/workflow/data/fasta/hs37d5.fasta

if [ "${BQSR,,}" = "true" ] ;
then

    echo "## START BASE RECALIBRATOR" 

    gatk BaseRecalibrator -I ${enddir}/${filename}.sort.marked_duplicates.bam -R ${REFERENCE_GENOME} --known-sites data/gnomad_v2.1_sv.sites.vcf.gz -O ${enddir}/${filename}.recal_data.table1 --tmp-dir ${TMP_DIR} --read-validation-stringency LENIENT --add-output-sam-program-record true 

    echo "## BASE RECALIBRATOR DONE" 

#####################
### ApplyBQSR
#####################

    echo "## START APPLY BQSR" 

    gatk ApplyBQSR -R ${REFERENCE_GENOME} -I ${enddir}/${filename}.sort.marked_duplicates.bam --bqsr-recal-file ${enddir}/${filename}.recal_data.table1 -O ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam  --add-output-sam-program-record true --create-output-bam-index false --read-validation-stringency LENIENT --tmp-dir ${TMP_DIR} 

#--create-output-bam-index true #--preserve-qscores-less-than 10 
    echo "## APPLY BQSR DONE" 

else
    echo "WARNING: BQSR is not applied" 
fi

###############################################################

#####################
## Remove intermediate files
#####################

echo "Removing intermediate files"

rm ${enddir}/${filename}.sort.bam

if [ "${RM_SORT_MARKED_FILE,,}" = "true" ] && [ "${BQSR,,}" = "true" ];
then

    rm ${enddir}/${filename}.sort.marked_duplicates.bam
    
fi

#rm ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam

################################################################

#######################
## Create Final Index
#######################

echo "Start create index" 
if [ "${RM_SORT_MARKED_FILE,,}" = "true" ] && [ "${BQSR,,}" = "true" ];
then
    samtools index ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam 
elif [ "${RM_SORT_MARKED_FILE,,}" = "false" ] && [ "${BQSR,,}" = "true" ];
then    
    samtools index ${enddir}/${filename}.sort.marked_duplicates.bam 
    samtools index ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam 
else
    samtools index ${enddir}/${filename}.sort.marked_duplicates*bam 
fi
echo "Index done"

###############################################################

# FINAL STATISTICS

if [ "${FINAL_STATISTICS,,}" = "true" ];
then
    
    echo "## FINAL STATISTICS"

    cd ${enddir}

    if [ "${BQSR,,}" = "true" ];
    then

        bash calculate_statistics.sh ${filename}.sort.marked_duplicates.BQSR 

        bash calculate_depth_coverage.sh ${filename}.sort.marked_duplicates.BQSR 

    else

        bash calculate_statistics.sh ${filename}.sort.marked_duplicates 

        bash calculate_depth_coverage.sh ${filename}.sort.marked_duplicates
    fi

    cd .. # Retrun to the original directory

    # How many records are in the filtered BAM compared to the original? How many read pairs does this represent?
    echo "How many records are in the filtered BAM compared to the original:" > ${enddir}/${filename}.records_last

    if [ "${BQSR,,}" = "true" ];
    then
        echo "`samtools view -c ${fullname}` `samtools view -c ${enddir}/${filename}.sort.marked_duplicates.BQSR.bam`" | awk '{printf "%0.2f%%\n", 100*$2/$1}' >> ${enddir}/${filename}.records_last
    else
        echo "`samtools view -c ${fullname}` `samtools view -c ${enddir}/${filename}.sort.marked_duplicates.bam`" | awk '{printf "%0.2f%%\n", 100*$2/$1}' >> ${enddir}/${filename}.records_last
    fi

    echo "## FINAL STATISTICS DONE"

fi

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

#echo "Start BaseRecalibrator after" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk BaseRecalibrator -I ${enddir}/${filename}.sort.marked_duplicates*HQ.bam -R ${REFERENCE_GENOME} --known-sites /vault/bio-scratch/arnau/structural_variants/gnomad_v2.1_sv.sites.vcf.gz -O ${enddir}/${filename}.recal_data.table2 && echo "End BaseRecalibrator after"

######################
### AnalyzeCovariates
#####################

#echo "Start AnalyzeCovariates" && /vault/mauricio/bio_team/gatk/gatk-4.2.6.0/./gatk AnalyzeCovariates -before ${enddir}/${filename}.recal_data.table1 -after ${enddir}/${filename}.recal_data.table2 -plots ${enddir}/${filename}.AnalyzeCovariates_last.pdf && echo "End AnalyzeCovariates"

##############################################################

# Remove temporary files

N_TMP_FILES=$(ls ${TMP_DIR} | wc -l) # Set number of temporary files

if [ "${RM_TMP_FILES,,}" = "true" ];
then
    if [[ $N_TMP_FILES -gt 0 ]];
    then
        echo "Removing all temporary files" &>> ${enddir}/${filename}.log
        rm -r ${TMP_DIR}/*
    fi
fi

# Change file permissions
#chmod ugo-x ${enddir}/${filename}.sort.marked_duplicates.bam

echo "${filename}.bam processing is finished"

done

###############################################################
###############################################################

echo "##### END OF THE PIPELINE #####" >&3
echo "Check the results in '${enddir}' folder!" >&3
