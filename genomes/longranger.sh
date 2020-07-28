#!/bin/bash

# Longranger variant calling script for 10x data

# Note: ensure to have followed the instructions on the website for creating the required folders and files and altering the reference fasta as needed

#Program versions tested with this script
#longranger v2.2.2

# Compute requirements:
# 256GB RAM
# 32 CPUs
# 175h Walltime
# 250GB Storage

# Programs that must be added to path
# longranger

#Set Required Variables
DIR=/path/to/output
OUTPUT=sample
REF=/path/to/refdatadir
FASTQS=/path/to/fastqs
SEX=female
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Increase file size limit
ulimit -n 8192

#Move into required directory
cd /mnt/volume1/antechinus_longranger

#Run long ranger WGS pipeline
longranger wgs \
--id=${SAMPLE} \
--reference=${REF} \
--fastqs=${FASTQS} \
--vcmode=freebayes \
--sex=${SEX}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "Long Ranger Pipeline complete. Please copy data to s3."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Long Ranger Pipeline complete. Please copy data to s3."}' ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error during Long Ranger pipeline. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during Long Ranger pipeline. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi