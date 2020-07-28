import os, sys, re, subprocess
exec(open('/usr/share/Modules/init/python.py').read())
module('load', 'bedtools')
from subprocess import Popen, PIPE

#This script takes a list of gene IDs and generates a fasta file needed for exon capture myBaits probe design
#written in python/2.7.9

#Usage: python generate_myBaits_targets_from_gene_ID.py genome_annotation_gtf list_of_gene_ID genome_assembly_fasta output_file_prefix

#Testing done on Artemis; will need to modify code to provide path to bedtools if run elsewhere

#Please contact Yuanyuan in case of any issues


gtf_file = sys.argv[1]
gtf_data = open(gtf_file).read().rstrip("\n").split("\n")

gene_ID_file = sys.argv[2]
gene_ID_list = open(gene_ID_file).read().rstrip("\n").split("\n")

genome_fasta_file = sys.argv[3]

outfile_prefix = sys.argv[4]
outfile_tmp_name = outfile_prefix + ".tmp"
outfile_tmp = open(outfile_tmp_name, "w")
outfile_bed_name = outfile_prefix + ".bed"
outfile_bed = open(outfile_bed_name, "w")
outfile_fasta_name = outfile_prefix + ".fasta"
outfile_summary_name = outfile_prefix + ".summary.txt"
outfile_summary = open(outfile_summary_name, "w")

gene_list = []
target_list = []

print("generating bed file...\n")

for line in gtf_data:
    if line.startswith("#"):
        continue

    line_split = line.split("\t")

    if line_split[2] != "exon":
        continue

    description = line_split[8]
    description_split = description.split(";")

    gene_id = []
    exon_id = []

    for subfield in description_split:
        if "gene_id" in subfield:
            gene_id = subfield.split("\"")[1]
        if "exon_id" in subfield:
            exon_id = subfield.split("\"")[1]

    if not gene_id:
        print "Warning: Missing gene_id in feature:", line
        #sys.exit()

    if gene_id in gene_ID_list:
        scaffold = line_split[0]
        target_start = int(line_split[3]) - 41
        target_end = int(line_split[4]) + 40
        strand = line_split[6]

        if not exon_id:
            target_id = gene_id
        else:
            target_id = "_".join([gene_id,exon_id])

        target = "\t".join([scaffold,str(target_start),str(target_end),target_id,"0",strand])

        if target not in target_list:
            target_list.append(target)
        if gene_id not in gene_list:
            gene_list.append(gene_id)

for target in target_list:
    outfile_tmp.write(target)
    outfile_tmp.write("\n")

outfile_tmp.close()

for gene in gene_ID_list:
    if gene not in gene_list:
        print "Warning: gene", gene, "not found in gtf file"

sort = ["sort", "-k1,1", "-k2,2n", outfile_tmp_name]
merge = ["bedtools", "merge", "-c", "4,5,6", "-o", "distinct", "-i", "-"]

proc1 = Popen(sort, stdout=PIPE)
proc2 = Popen(merge, stdin=proc1.stdout, stdout=outfile_bed)
proc2.wait()

outfile_bed.close()
os.remove(outfile_tmp_name)

print("\nextracting target sequences...\n")

getfasta = ["bedtools", "getfasta", "-s", "-fi", genome_fasta_file, "-bed", outfile_bed_name, "-fo", outfile_fasta_name]

proc2 = Popen(getfasta)
proc2.wait()
if proc2.returncode < 0:
    print "Error: bedtools crashed with signal", -proc.returncode, "(Note: This script can't run on HPC login node. Use interactive or batch job)"
    sys.exit()

print("\nwriting summary file...\n")

def line_count(file):
    count = 0
    for line in open(file):
        count += 1
    return count

def seq_count(fasta_file):
    count = 0
    for line in open(fasta_file):
        if line.startswith(">"):
            count += 1
    return count

outfile_summary.write("Number of targeted genes: %d%s" % (len(gene_list),"\n"))
outfile_summary.write("Number of targets: %d%s" % (line_count(outfile_bed_name),"\n"))
outfile_summary.write("Number of target sequences: %d%s" % (seq_count(outfile_fasta_name),"\n"))

outfile_summary.close()

print("\nrun completed")
