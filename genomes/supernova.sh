#!/bin/bash

# Supernova de novo genome assembly script

#Program versions tested with this script
#supernova v2.1.1(6bb16452a)

# Compute requirements:
# 768GB RAM
# 96 CPUs
# 150h Walltime
# 2.5TB Storage

# Programs that must be added to path
# supernova

#Set Required Variables
TOP_DIR="/path/to/dir/containing/fastqfolder"
FASTQ_DIR="${OUTPUT_DIR}/fastqfolder"
RUN_ID="species"
CPUS="96"
MEM="768"
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

# Move into required directory which contains folder with input fastqs (DATASET_NAME) and where all outputs will be written
cd ${TOP_DIR}

# Run the Supernova assembly
date
echo "Starting supernova de novo genome assembly..."
supernova run \
    --id=${RUN_ID} \
    --fastqs=${FASTQ_DIR} \
    --maxreads=all \
    --localcores=${CPUS} \
    --localmem=${MEM}

# If supernova ran without errors, move onto creating output FASTAs, otherwise terminate VM
STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "Supernova run complete. Creating output FASTA files..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Supernova run complete. Creating output FASTA files..."}' ${SLACK_URL}
    else
        date
        echo "Error during supernova run. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during supernova run. Please check VM. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

# Generate FASTA outputs
date
echo "Starting creation of supernova FASTA output files..."
supernova mkoutput \
    --style=pseudohap2 \
    --asmdir=${RUN_ID}/outs/assembly \
    --outprefix=${RUN_ID}_pseudohap2 \
    --index

supernova mkoutput \
    --style=raw \
    --asmdir=${RUN_ID}/outs/assembly \
    --outprefix=${RUN_ID}_raw

supernova mkoutput \
    --style=megabubbles \
    --asmdir=${RUN_ID}/outs/assembly \
    --outprefix=${RUN_ID}_megabubbles

supernova mkoutput \
    --style=pseudohap \
    --asmdir=${RUN_ID}/outs/assembly \
    --outprefix=${RUN_ID}_pseudohap

# Check to see if final files have been created and terminate the VM
if [[ -s ${RUN_ID}_pseudohap.fasta.gz ]];
    then
        date
        echo "Final FASTA files created. Please copy data to a safe place. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Final FASTA files created. Please copy data to a safe place. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error when creating final FASTA files. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error when creating final FASTA files. Please check VM. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi