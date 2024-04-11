#!/usr/bin/python



import os, re, sys, string

# Use 'find . -maxdepth 3 -name "fastq" > read_files' to identify the folders and read files of interest
file_list = open("read_files", "r")
#ref = sys.argv[1]
#bowtieref = ref.replace(".fasta","")

# Counter to figure out the number of jobs run
count = 0

for line in file_list:
    #pull out the line and strip off the newline and extra stuff
    path = line.strip()
    directory = path.replace(".fastq","").split("/")[0]
    name = path.replace(".fastq","").split("/")[1]
    name1 = path.replace(".fastq","").split("/")[1]
    spec = path.split(".")[1].split("/")[0]
    print(directory)
    count += 1
    print("Processing" + name + " now...")

    if os.path.isfile(name + "_job.sh")==True:
        name = name + "_lym"

    # GENERATE COMMANDS FOR BASH SCRIPTS
    ######### - FILTER AND TRIM - #########
    # Start quality filter
    qual = "QualFilterFastq.pl -i " + name1 + ".fastq -m 20 -x 10 -o " + name + "_trimmed.fastq"
    # Homopolymer repeat filter
    HRF = "HRFilterFastq.pl -i " + name + "_trimmed.fastq -n 30 -o " + name + "_trimmedfiltered.fastq"
    # Adaptor filter
    adapter = "bbduk.sh in=" + name + "_trimmedfiltered.fastq ref=/home/lmesrop/miniconda3/opt/bbmap-38.18/resources/adapters.fa k=12 stats=" + name + "stats.txt out=" + name + "_trimfiltclean.fastq overwrite=t"
    # Remove PCR duplicates
    dups = "RemovePCRDups.pl -i " + name + "_trimfiltclean.fastq -o " + name + "_trimfiltcleanrem.fastq -s 1 -e 4 -j 9"
    # Trim tags
    tags = "TagTrimmer.pl -i " + name + "_trimfiltcleanrem.fastq -b 1 -e 8 -o " + name + "_trimfiltcomplete.fastq"

    ######### - READ MAPPING - #########
    # Align reads to reference
    align = "bowtie2 -x " + "/home/lmesrop/ref_transcriptomes/Vargula_tsujii_cdhit_95_bowtie_index/Vargula_tsujii_cdhit_95" + " -U " + name + "_trimfiltcomplete.fastq -S " + name + "_bowtie2.sam"
    # Filter mapped reads to remove unmapped and non-primary
    second = "samtools view -hb -F 2308 " + name + "_bowtie2.sam | samtools sort - > " + name + "_bowtie2.sorted.bam"
    # Re-format bam file to sam
    reformat = "samtools view -h -o " + name + "_bowtie2.sorted.sam " + name + "_bowtie2.sorted.bam"

    # Filter alignments for weak/ambiguous/duplicated alignments
    filt = "SAMFilterByGene.pl -i " + name + "_bowtie2.sorted.sam -m 40 -o " + name + "_final.sam -p 1"

    os.system("echo \"#!/bin/bash\n\n#PBS -N taqseq\n" +  "#PBS -l nodes=1:ppn=4\n#PBS -l walltime=108:00:00\n#PBS -m e -M lmesrop@ucsb.edu\n\n" + "cd $PBS_O_WORKDIR\n\n"  + "source /home/lmesrop/.bash_profile\n\n" + qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align +
"\n" + second + "\n" + reformat + "\n" + filt + "\" > " + name + "_job.sh")

    os.system("pbs" + name + "_job.sh")

file_list.close()
print("Completed processing" + str(count) + " files.")

#qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align + "\n" + filt + "\" > " + name + "_job.sh")


