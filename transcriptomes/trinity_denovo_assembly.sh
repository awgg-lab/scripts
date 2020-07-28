#!/bin/bash

# Trinity de novo transcriptome assembly script (including trimmomatic trimming of reads)

#Program versions tested with this script
#trinity v2.8.6 (all correct versions of dependencies will be installed with docker)

# Compute requirements:
# 384GB RAM
# 48 CPUs
# 100h Walltime
# 750GB Storage

# This script must be run as the root user to enable trinity to be installed via docker.
# Otherwise hash out the docker and trinity installation steps below

#Install Docker
yes | apt update
yes | apt install docker.io
yes | systemctl start docker
yes | systemctl enable docker

#Pull Trinity Container
yes | docker pull trinityrnaseq/trinityrnaseq

# Set Variables
FASTQ_DIR="/path/to/fastqs"
OUTPUT="species_sample_trinity"
R1="sample_R1_trimmed_P.fastq.gz"
R2="sample_R2_trimmed_P.fastq.gz"
#SAMPLES="trinity_samples_file.txt" #Use this option instead of R1 and R2 if you have many samples/read files
CPUS="48"
MEM="384"
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Run Trinity pipeline
date
echo "Starting Trinity de novo transcriptome assembly..."
docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
--trimmomatic \
--seqType fq --max_memory ${MEM} \
--left `pwd`/${R1} \
--right `pwd`/${R2} \
#--samples_file `pwd`/${SAMPLES} #Use this option instead of --left and --right if you have many samples/read files
--output `pwd`/${output} \
--CPU ${CPUS} \
--full_cleanup \

#Check if Trinity completed successfully, if so continue to calculate basic summary stats, otherwise terminate the VM.
if [[ -s ${OUTPUT}.Trinity.fasta ]];
    then
        date
        echo "Trinity assembly complete. Running basic assembly stats..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Trinity assembly complete. Running basic assembly stats..."}' \
        ${SLACK_URL}
    else
        date
        echo "Error when running trinity assembly, please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running trinity assembly, please check VM. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Run basic stats on output assembly
docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq \
/usr/local/bin/trinityrnaseq/util/TrinityStats.pl \
`pwd`/${OUTPUT}.Trinity.fasta > `pwd`/${OUTPUT}.Trinity.stats

#Check that the assembly stats have been saved and terminate the VM.
if [[ -s ${OUTPUT}.Trinity.stats ]];
    then
        date
        echo "Trinity pipeline complete. Please copy data to a safe place. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Trinity pipeline complete. Please copy data to a safe place. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error when running trinity stats, please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running trinity stats, please check VM. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi