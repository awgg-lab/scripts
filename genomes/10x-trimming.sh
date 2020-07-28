#!/bin/bash

# 10X Trimming script to convert 10x genomics reads to standard illumina reads

#Program versions tested with this script
#bbmap v38.73

# Compute requirements:
# 512GB RAM
# 64 CPUs
# 5h Walltime
# 500GB Storage (~3x the size of the input reads)

# Programs that must be added to path
# bbduk.sh
# repair.sh

#Set Required Variables
FQPATH=/path/to/fastqs
OUTPUT=sample
CPUS=64
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Move into input directory
cd $FQPATH

#Trim first 23 bases off R1 Reads - 16bp barcode + 7bp poor sequence quality
bbduk.sh -Xmx432g in=${OUTPUT}_R1.fastq.gz out=${OUTPUT}_R1.trim.fq.gz ftl=23 ordered

if [[ -s ${OUPUT}_R1.trim.fq.gz ]];
    then
        date
        echo "R1 trimming complete. Moving on to R2 trimming..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"R1 trimming complete. Moving on to R2 trimming..."}' ${SLACK_URL}
    else
        date
        echo "Error during R1 trimming. Please check VM."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during R1 trimming. Please check VM."}' ${SLACK_URL}
fi

#Trim first base off R2 reads due to poor sequence quality
bbduk.sh -Xmx432g in=${OUTPUT}_R2.fastq.gz out=${OUTPUT}_R2.trim.fq.gz ftl=1 ordered

if [[ -s ${OUTPUT}_R2.trim.fq.gz ]];
    then
        date
        echo "R2 trimming complete. Checking pairing..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"R2 trimming complete. Checking pairing..."}' ${SLACK_URL}
    else
        date
        echo "Error during R2 trimming. Please check VM."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during R2 trimming. Please check VM."}' ${SLACK_URL}
fi

#Check paired end reads match up correctly after trimming
repair.sh t=$CPUS in1=${OUTPUT}_R1.trim.fq.gz in2=${OUTPUT}_R2.trim.fq.gz out1=${OUTPUT}_R1.trim.pair.fq.gz out2=${OUTPUT}_R2.trim.pair.fq.gz

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "10X Trimming pipeline complete. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"10X Trimming pipeline complete. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error when checking pairing. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error when checking pairing. Please check VM. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi