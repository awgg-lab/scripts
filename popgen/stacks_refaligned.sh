#!/bin/bash

# Stacks refaligned script

# Compute requirements:
# 256GB RAM
# 32 CPUs
# 10h Walltime
# ~1GB storage per sample

# Program versions for this script
# aws-cli=1.16.152
# bwa=0.7.17-r1188
# samtools=1.7
# stacks=2.53
# GNU parallel=20161222

# User Defined Variables
TOPDIR="/mnt/volume1"
LANE1="${TOPDIR}/raw/lane"
BARCODES1="${TOPDIR}/raw/lane/barcodes.txt"
GENOME="${TOPDIR}/genome/genome.fa"
CLEANED="${TOPDIR}/cleaned"
REFALIGNED="${TOPDIR}/refaligned"
OUTDIR="${TOPDIR}/gstacks"
POPMAP="${TOPDIR}/gstacks/popmap.txt"
MIN_SAMPLES="0.3"
MIN_MAF="0.01"
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

# Automatically set the number of CPUs
CPUS=$(nproc --all)

# Change the limit on the number of files that can be open at once (in case you have a large number of samples)
ulimit -n 8192

#Move into top directory
cd ${TOPDIR}

#Step 1: Remove barcodes
date
echo "Starting Process radtags..."
process_radtags -c -q -i gzfastq -p ${LANE1} --disable_rad_check --inline_null -b ${BARCODES1} -o ${CLEANED}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]]; then
    date
    echo "Process radtags complete. Preparing reference genome..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Process radtags complete. Preparing reference genome..."}' \
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

#Step 2: Check if genome is indexed and if not, index it
date
if [[ -s ${GENOME}.bwt ]]; then
    date
    echo "Genome already indexed. Aligning reads..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Genome already indexed Aligning reads..."}' \
    ${SLACK_URL}
else
    echo "Starting genome indexing..."
    bwa index ${GENOME}
    STATUS=$(echo $?)
    if [[ ${STATUS} -eq 0 ]]; then
        date
        echo "Genome indexing complete. Aligning reads..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Genome indexing complete. Aligning reads..."}' \
        ${SLACK_URL}
    else
        date
        echo "Error during genome indexing. Please check VM. Stopping instance..."
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error during genome indexing. Please check VM. Stopping instance..."}' \
        ${SLACK_URL}
        exit 1
    fi
fi

#Step 3: Align raw reads to reference genome. Convert .sai to .bam files and sort .bam files for gstacks input
date
echo "Starting alignment of reads to the reference genome..."

for i in $(cut -f 1 ${POPMAP}); do
    if [[ -s ${TOPDIR}/refaligned/${i}.sort.bam ]]; then
        echo "Sample ${i} already aligned. Skipping..."
    else
        echo "${i}" >> ${TOPDIR}/align_samples.txt
    fi
done

ALIGN_SAMPLES="${TOPDIR}/align_samples.txt"

for sample in $(cat ${ALIGN_SAMPLES}); do
    bwa aln -B 5 -t ${CPUS} ${GENOME} ${CLEANED}/${sample}.fq.gz > $REFALIGNED/${sample}.sai
done

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]]; then
    date
    echo "Reads successfully aligned to reference. Converting and sorting output..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Reads successfully aligned to reference. Converting and sorting output..."}' \
    ${SLACK_URL}
else
    date
    echo "Error when mapping reads to reference. Please check VM. Stopping instance..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Error when mapping reads to reference. Please check VM. Stopping instance..."}' \
    ${SLACK_URL}
    shutdown -h 15
    exit 1
fi

date
echo "Starting conversion and sorting of alignments..."

parallel "bwa samse ${GENOME} ${REFALIGNED}/{}.sai ${CLEANED}/{}.fq.gz | samtools view -b > ${REFALIGNED}/{}.bam" < ${ALIGN_SAMPLES}

for sample in $(cat ${ALIGN_SAMPLES}); do
    samtools sort -@ ${CPUS} ${REFALIGNED}/${sample}.bam -o ${REFALIGNED}/${sample}.sort.bam
done

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]]; then
    date
    echo "Ouput converted and sorted successfully. Running gstacks..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Ouput converted and sorted successfully. Running gstacks..."}' \
    ${SLACK_URL}
else
    date
    echo "Error when converting and sorting output. Please check VM. Stopping instance..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Error when converting and sorting output. Please check VM. Stopping instance..."}' \
    ${SLACK_URL}
    shutdown -h 15
    exit 1
fi

#Step 4: Run gstacks pipeline to call SNPs
date
echo "Starting gstacks..."
gstacks -I ${REFALIGNED} -O ${OUTDIR} -M ${POPMAP} -S .sort.bam -t ${CPUS}

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]]; then
    date
    echo "Gstacks complete. Running populations module..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Gstacks complete. Running populations module..."}' \
    ${SLACK_URL}
else
    date
    echo "Error during gstacks. Please check VM. Stopping instance..."
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Error during gstacks. Please check VM. Stopping instance..."}' \
    ${SLACK_URL}
    shutdown -h 15
    exit 1
fi

#Step 5: Run populations module to filter and export genotype calls
date
echo "Starting populations..."
populations -P ${OUTDIR} -M ${POPMAP} -r ${MIN_SAMPLES} --min-maf ${MIN_MAF} --write-random-snp --vcf -t ${CPUS} &> ${OUTDIR}/populations.vcf.oe

STATUS=$(echo $?)
if [[ ${STATUS} -eq 0 ]]; then
    echo "Stacks pipeline complete. Stopping instance..."
    date
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Stacks pipeline complete. Stopping instance..."}' \
    ${SLACK_URL}
    shutdown -h 15
    exit 0
else
    echo "Error during population module. Please check VM. Stopping instance..."
    date
    curl -X POST -H 'Content-type: application/json' \
    --data '{"text":"Error during population module. Please check VM. Stopping instance..."}' \
    ${SLACK_URL}
    shutdown -h 15
    exit 1
fi