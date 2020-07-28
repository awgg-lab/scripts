#!/bin/bash

#Script to run RepeatModeler to create a custom repeat database for genome annotation

#Program versions tested with this script
#repeatmodeler v2.0.1

#Compute requirements
# 64GB RAM
# 16 CPUs
# 15h Walltime
# 50GB Storage

#Set variables
DIR=/path/to/input
REFERENCE=ref.fasta
OUTPUT=output_prefix
CPUS=16
PA=`expr $CPUS - 2`
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Move into maker directory
cd ${DIR}

#Begin script
echo "RepeatModeler log file"
pwd
date

#Build repeat database
BuildDatabase -name ${OUTPUT}_repeat_db -engine ncbi ${REFERENCE}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "Repeat database built. Starting RepeatModeler..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Repeat database built. Starting RepeatModeler..."}' \
        ${SLACKURL}
    else
        date
        echo "Error when building repeat database, please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when building repeat database, please check VM. Stopping instance..."}' \
        ${SLACKURL}
        shutdown -h 15
        exit 1
fi

#Run RepeatModeler
RepeatModeler -pa ${PA} -engine ncbi -database ${OUTPUT}_repeat_db

STATUS2=$(echo $?)
if [[ ${STATUS2} -eq 0 ]];
    then
        date
        echo "RepeatModeler pipeline complete. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"RepeatModeler pipeline complete. Stopping instance..."}' \
        ${SLACKURL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error during RepeatModeler pipeline. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error during RepeatModeler pipeline. Please check VM. Stopping instance..."}' \
        ${SLACKURL}
        shutdown -h 15
        exit 1
fi