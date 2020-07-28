#!/bin/bash

#Script to run RepeatMasker to mask a genome using a custom repeat database from RepeatModeler 

#Program versions tested with this script
#repeatmasker v4.0.6

#Compute requirements
# 64GB RAM
# 16 CPUs
# 15h Walltime
# 10GB Storage

#Set variables
DIR=/path/to/output
REF=/path/to/ref.fasta
REPEATLIB=/path/to/consensi.fa.classified
CPUS=16
PA=`expr $CPUS - 2`
SLACKURL="${SLACK_URL}"SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Begin script
echo "RepeatMasker log file"
pwd
date

#Run RepeatMasker
RepeatMasker -dir ${DIR} -pa ${PA} -gff -lib ${REPEATLIB} -nolow ${REF}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "RepeatMasker Complete. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"RepeatMasker Complete. Stopping instance..."}' ${SLACKURL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error during RepeatMasker. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during RepeatMasker. Please check VM. Stopping instance..."}' ${SLACKURL}
        shutdown -h 15
        exit 1
fi