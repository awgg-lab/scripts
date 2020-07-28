#!/bin/bash

#Trinotate transcriptome annotation pipeline

#Program versions tested with this script
#Trinotate v3.2.0
#TransDecoder v5.5.0
#RNAMMER v1.2
#HMMER v3.2.0
#SignalP v4.1
#SQLite v3.29

# Compute requirements:
# 128GB RAM
# 16 CPUs
# 150h Walltime
# 100GB Storage

#Ensure trinotate and all dependencies are installed and in the current path
#Ensure trinotate_db containing required databses has been made available - see trinotate documentation for more information

#Set Variables
TRANSCRIPTOME="/path/to/sample.Trinity.fasta"
DATABASE="/path/to/trinotate_db"
OUT_DIR="path/to/outputdir"
RNAMMER="path/to/rnammer"
CPUS="16"
MEM="128"
SLACK_URL="https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX"

#Move into output directory
cd ${OUT_DIR}

#Get sample name from trinity fasta file
SAMPLE=`basename -s .Trinity.fasta ${TRANSCRIPTOME}`

#Copy Trinotate.sqlite file from database to current folder
cp ${DATABASE}/Trinotate.sqlite .

#Step 1: TransDecoder
echo "Starting Step 1: TransDecoder $(date)"

#extract the long open reading frames, deafult is ORFs >100aa long
TransDecoder.LongOrfs -t ${TRANSCRIPTOME}

#predict likely coding regions
TransDecoder.Predict -t ${TRANSCRIPTOME}

if [[ -s ${TRANSCRIPTOME}.transdecoder.pep ]];
    then
        echo "TransDecoder step complete. Moving on to RNAMMER..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"TransDecoder step complete. Moving on to RNAMMER..."}' \
        ${SLACK_URL}
    else
        echo "Error when running TransDecoder step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running TransDecoder step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 2: RNAMMER
echo "Starting Step 2: RNAMMER $(date)"

#Concatenate all transcripts together into a super-scaffold, run RNAMMER to ID rRNA homologies, then transform the rRNA coordinates in the superscaffold back to the transcriptome reference coordinates.
perl /apps/Trinotate-v3.2.0/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ${TRANSCRIPTOME} --path_to_rnammer ${RNAMMER}

if [[ -s ${TRANSCRIPTOME}.rnammer.gff ]];
    then
        echo "RNAMMER step complete. Moving on to BLAST..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"RNAMMER step complete. Moving on to BLAST..."}' \
        ${SLACK_URL}
    else
        echo "Error when running RNAMMER step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running RNAMMER step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 3: BLAST
echo "Starting Step 3: BLAST $(date)"

#1. search trinity transcripts against swissprot (don't change output file name after >)
blastx -query ${TRANSCRIPTOME} -db ${DATABASE}/uniprot_sprot.pep -num_threads ${CPUS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastx_1e-5.outfmt6

#2. search transdecoder predicted proteins against swissprot (don't change output file name after >)
blastp -query ${TRANSCRIPTOME}.transdecoder.pep -db ${DATABASE}/uniprot_sprot.pep -num_threads ${CPUS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_1e-5.outfmt6

#3. search transdecoder predicted proteins against devil genome annotations (don't change output file name after >)
blastp -query ${TRANSCRIPTOME}.transdecoder.pep -db ${DATABASE}/Sarcophilus_harrisii.DEVIL7.0.pep.all.fa -num_threads ${CPUS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_devil_1e-5.outfmt6

#4. search trinity transcripts against IDMM (don't change output file name after >)
blastx -query ${TRANSCRIPTOME} -db ${DATABASE}/allSeq.fasta -num_threads ${CPUS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastx_1e-5_IDMM.outfmt6

#5. search transdecoder predicted proteins against IDMM (don't change output file name after >)
blastp -query ${TRANSCRIPTOME}.transdecoder.pep -db ${DATABASE}/allSeq.fasta -num_threads ${CPUS} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_1e-5_IDMM.outfmt6

COUNT=$(ls -l *.outfmt6 | wc -l)
if [[ ${COUNT} -eq 5 ]];
    then
        echo "BLAST step complete. Moving on to hmmscan..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"BLAST step complete. Moving on to hmmscan..."}' \
        ${SLACK_URL}
    else
        echo "Error when running BLAST step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running BLAST step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 4: hmmscan
echo "Starting Step 4: hmmscan $(date)"

#run hmmscan to identify protein families within predicted protein sequences. This script assumes you've already run transdecoder. Don't change the name of the output file (after >)
hmmscan --cpu ${CPUS} --domtblout TrinotatePFAM.out ${DATABASE}/Pfam-A.hmm ${TRANSCRIPTOME}.transdecoder.pep > pfam.log

if [[ -s TrinotatePFAM.out ]];
    then
        echo "hmmscan step complete. Moving on to signalP..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"hmmscan step complete. Moving on to signalP..."}' \
        ${SLACK_URL}
    else
        echo "Error when running hmmscan step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running hmmscan step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 5: signalP
echo "Starting Step 5: signalP $(date)"

#run signalP
signalp -f short -n signalp.out ${TRANSCRIPTOME}.transdecoder.pep

if [[ -s signalp.out ]];
    then
        echo "signalP step complete. Moving on to Trinotate..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"signalP step complete. Moving on to Trinotate..."}' \
        ${SLACK_URL}
    else
        echo "Error when running signalP step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running signalP step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 6: Trinotate
echo "Starting Step 6: Trinotate $(date)"

#import the transcriptome (Trinity.fasta), protein sequences (transdecoder) and gene/transcript map file (above) into SQLite
Trinotate Trinotate.sqlite init \
--gene_trans_map ${TRANSCRIPTOME}.gene_trans_map \
--transcript_fasta ${TRANSCRIPTOME} \
--transdecoder_pep ${TRANSCRIPTOME}.transdecoder.pep

#load BLAST search results
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp_1e-5.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx_1e-5.outfmt6
Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 blastp_devil_1e-5.outfmt6 --prog blastp --dbtype ${DATABASE}/Sarcophilus_harrisii.DEVIL7.0.pep.all.fa
Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 blastp_1e-5_IDMM.outfmt6 --prog blastp --dbtype ${DATABASE}/allSeq.fasta
Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 blastx_1e-5_IDMM.outfmt6 --prog blastx --dbtype ${DATABASE}/allSeq.fasta

#load pfam domain results
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

# load signal peptide predictions
Trinotate Trinotate.sqlite LOAD_signalp signalp.out

#load rnammer results
Trinotate Trinotate.sqlite LOAD_rnammer ${TRANSCRIPTOME}.rnammer.gff

#generate trinotate report - change the .xls file name as appropriate
Trinotate Trinotate.sqlite report -E 1e-5 --pfam_cutoff DNC --incl_pep --incl_trans > ${SAMPLE}_annotation_report.xls

if [[ -s ${SAMPLE}_annotation_report.xls ]];
    then
        echo "Trinotate step complete. Moving on to extract GO terms..."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Trinotate step complete. Moving on to extract GO terms..."}' \
        ${SLACK_URL}
    else
        echo "Error when running Trinotate step. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when running Trinotate step. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi

#Step 7: extract GO terms
echo "Starting Step 7: Extract GO Terms $(date)"

#extract GO terms from trinotate.xls file
perl /apps/Trinotate-v3.2.0/util/extract_GO_assignments_from_Trinotate_xls.pl \
--Trinotate_xls ${SAMPLE}_annotation_report.xls -T -I > ${SAMPLE}_go_inclancestral.txt

if [[ -s ${SAMPLE}_go_inclancestral.txt ]];
    then
        echo "Trinotate pipeline complete! :-) Please transfer data to s3"
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Trinotate pipeline complete! :-) Please transfer data to s3"}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 0
    else
        echo "Error when extracting GO terms. Please check log files."
        date
        curl -X POST -H 'Content-type: application/json' \
        --data '{"text":"Error when extracting GO terms. Please check log files."}' \
        ${SLACK_URL}
        shutdown -h 15
        exit 1
fi