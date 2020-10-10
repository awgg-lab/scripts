#!/bin/bash

# Stacks de novo script

# Compute requirements:
# 256GB RAM
# 32 CPUs
# 60h Walltime
# ~1GB storage per sample

# User Defined Variables
TOPDIR="/mnt/volume1"
LANE1="${TOPDIR}/raw/lane"
BARCODES1="${TOPDIR}/raw/lane/barcodes.txt"
CLEANED="${TOPDIR}/cleaned"
DENOVO="${TOPDIR}/denovo"
OUTDIR="${TOPDIR}/gstacks"
POPMAP="${TOPDIR}/gstacks/popmap.txt"
MAX_DISTANCE="2"
MIN_DEPTH="3"
MISMATCHES="2"
MIN_SAMPLES="0.3"
MIN_MAF="0.01"
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"
SAMPLES=$(cut -f 1 ${POPMAP})

# Automatically set the number of CPUs
CPUS=$(nproc --all)

#Change the limit on the number of files that can be open at once due to large number of samples
ulimit -n 8192

#Move into top directory
cd ${TOPDIR}

#Step 1: Run process_radtags
date
echo "Starting Process radtags..."
process_radtags -c -q -i gzfastq -p ${LANE1} --disable_rad_check --inline_null -b ${BARCODES1} -o ${CLEANED}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "Process radtags complete. Starting ustacks..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Process radtags complete. Starting ustacks..."}' \
        ${SLACK_URL}
    else
        date
        echo "Error during process radtags. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error during process radtags. Please check VM. Stopping instance..."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 2: Run ustacks
date
echo "Starting ustacks ..."

id=1
for sample in ${SAMPLES}; do
    ustacks -f ${CLEANED}/${sample}.fq.gz -o ${DENOVO} -i $id -M ${MAX_DISTANCE} -m ${MIN_DEPTH} -p ${CPUS}
    let "id+=1"
done

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "ustacks complete. Moving on to cstacks..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"ustacks complete. Moving on to cstacks..."}' ${SLACK_URL}
    else
        date
        echo "Error while running ustacks. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error while running ustacks. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 3: Run cstacks
date
echo "Starting cstacks ..."
cstacks -n ${MISMATCHES} -P ${DENOVO} -M ${POPMAP} -p ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "cstacks complete. Moving on to sstacks..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"cstacks complete. Moving on to sstacks..."}' ${SLACK_URL}
    else
        date
        echo "Error during cstacks. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during cstacks. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 4: Run sstacks
date
echo "Starting sstacks ..."
sstacks -P ${DENOVO} -M ${POPMAP} -p ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "sstacks complete. Moving on to tsv2bam..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"sstacks complete. Moving on to tsv2bam..."}' ${SLACK_URL}
    else
        date
        echo "Error during sstacks. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during sstacks. Stopping instance...   "}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 5: Run tsv2bam to transpose the data so it is stored by locus, instead of by sample.
date
echo "Transposing data..."
tsv2bam -P ${DENOVO} -M ${POPMAP} -t ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "tsv2bam complete. Moving on to gstacks..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"tsv2bam complete. Moving on to gstacks..."}' ${SLACK_URL}
    else
        date
        echo "Error during tsv2bam. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during tsv2bam. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 4: Run gstacks
date
echo "Starting gstacks ..."
gstacks -P ${DENOVO} -O ${OUTDIR} -M ${POPMAP} -t ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "gstacks complete. Moving on to populations..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"gstacks complete. Moving on to populations..."}' ${SLACK_URL}
    else
        date
        echo "Error during gstacks. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during gstacks. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 4: Run populations
date
echo "Starting populations ..."
populations -P ${OUTDIR} -M ${POPMAP} -r ${MIN_SAMPLES} --write-single-snp --plink --vcf -t ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]];
    then
        date
        echo "Stacks pipeline complete. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Stacks pipeline complete. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        date
        echo "Error during populations. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' --data '{"text":"Error during populations. Stopping instance..."}' ${SLACK_URL}
        shutdown -h 15
        exit 1
fi