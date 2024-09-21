#!/usr/bin/python

#Description: Generate slurm submission scripts to process and map 3'Tag-Seq samples vtu1.counts.tab, vtu2.counts.tab, vtu3.counts.tab, vtu4.counts.tab, vtu5.counts.tab, vtu6.counts.tab, vtu7.counts.tab, vtu8.counts.tab, vtu9.counts.tab, vtu10.counts.tab. 
#Author: Jessica Goodheart 
#Goal: Take the raw fastq files, perform QC checks, filter the adaptor, remove PCR duplicates, trim, and map reads to reference transcriptome. 

import os, re, sys, string

# Use 'find . -maxdepth 3 -name "fastq" > read_files' to identify the folders and read files of interest
file_list = open("read_files", "r")
ref = sys.argv[1]
bowtieref = ref.replace(".fasta","")

# Counter to figure out the number of jobs run
count = 0

for line in file_list:
    #pull out the line and strip off the newline and extra stuff
    path = line.strip()
    directory = path.replace(".fastq","").split("/")[0]
    name = path.replace(".fastq","").split("/")[1]
    name1 = path.replace(".fastq","").split("/")[1]
    spec = path.split(".")[1].split("/")[0]
    print directory
    count += 1
    print "Processing " + name + " now..."

    if os.path.isfile(name + "_job.sh")==True:
        name = name + "_cdhit95"

    # GENERATE COMMANDS FOR BASH SCRIPTS
    ######### - FILTER AND TRIM - #########
    # Start quality filter
    qual = "QualFilterFastq.pl -i " + name1 + ".fastq -m 20 -x 10 -o " + name + "_trimmed.fastq"
    # Homopolymer repeat filter
    HRF = "HRFilterFastq.pl -i " + name + "_trimmed.fastq -n 30 -o " + name + "_trimmedfiltered.fastq"
    # Adaptor filter
    adapter = "bbduk.sh in=" + name + "_trimmedfiltered.fastq ref=../adapters.fasta k=12 stats=" + name + "stats.txt out=" + name + "_trimfiltclean.fastq overwrite=t"
    # Remove PCR duplicates
    dups = "RemovePCRDups.pl -i " + name + "_trimfiltclean.fastq -o " + name + "_trimfiltcleanrem.fastq -s 1 -e 4 -j 9"
    # Trim tags
    tags = "TagTrimmer.pl -i " + name + "_trimfiltcleanrem.fastq -b 1 -e 8 -o " + name + "_trimfiltcomplete.fastq"

    ######### - READ MAPPING - #########
    # Align reads to reference
    align = "bowtie2 -x " + bowtieref + " -U " + name + "_trimfiltcomplete.fastq -S " + name + "_bowtie2.sam"
    # Filter mapped reads to remove unmapped and non-primary
    second = "samtools view -hb -F 2308 " + name + "_bowtie2.sam | samtools sort - > " + name + "_bowtie2.sorted.bam"
    # Re-format bam file to sam
    reformat = "samtools view -h -o " + name + "_bowtie2.sorted.sam " + name + "_bowtie2.sorted.bam"

    # Filter alignments for weak/ambiguous/duplicated alignments
    filt = "SAMFilterByGene.pl -i " + name + "_bowtie2.sorted.sam -m 40 -o " + name + "_final.sam -p 1"

    os.system("echo \"#!/bin/bash\n\n#SBATCH --job-name=tagseq\n#SBATCH --workdir=/home/goodheart/projects/upper_lip_expression/tag-seq_data/v.tsujii/" + directory + "/ \n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=64G\n#SBATCH --time=72:00:00 \n\n" + qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align + "\n" + second + "\n" + reformat + "\n" + filt + "\" > " + name + "_job.sh")

    os.system("sbatch " + name + "_job.sh")

file_list.close()
print "Completed processing " + str(count) + " files."

#qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align + "\n" + filt + "\" > " + name + "_job.sh")
