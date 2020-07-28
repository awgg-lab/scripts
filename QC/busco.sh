#!/bin/bash

#Script to install and run BUSCO

#Program versions tested with this script
#busco v4.0.5 (all correct versions of dependencies will be isntalled with docker)

#Compute requirements
# 32GB RAM
# 8 CPUs
# 50h Walltime
# 15GB Storage

# This script must be run as the root user to enable BUSCO to be installed via docker.
# Otherwise hash out the installation steps below

#Install Docker
yes | apt update
yes | apt install docker.io
yes | systemctl start docker
yes | systemctl enable docker

#Pull BUSCO Container
yes | docker pull ezlabgva/busco:v4.0.5_cv1

#Assign variables
DIR=/path/to/inputdir
LINEAGE=mammalia_odb10 #run busco --list-datasets to see a complete list of available options
MODE=genome #Can also be transcriptome or proteins depending on input file
INPUT=file.fasta
OUTPUT=outputname
CPUS=8
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Run BUSCO
cd ${DIR}
date
echo "Starting BUSCO"

docker run --rm -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.5_cv1 busco -m ${MODE} -i ${INPUT} -o ${OUTPUT} -l ${LINEAGE} -c ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "BUSCO complete. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"BUSCO complete. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error during BUSCO. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error during BUSCO. Please check VM. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi