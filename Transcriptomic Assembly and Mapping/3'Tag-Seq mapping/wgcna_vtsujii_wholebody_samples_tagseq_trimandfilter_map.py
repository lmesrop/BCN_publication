#Update

#!/usr/bin/python

#Description: Generate slurm submission scripts to process and map 3'Tag-Seq samples A.1.counts.tab, A.2.counts.tab, A.3.counts.tab, AIF.1.counts.tab, AIF.2.counts.tab, AIF.3.counts.tab, AII.1.counts.tab, AII.2.counts.tab, AII.3.counts.tab, AIII.1.counts.tab, AIII.2.counts.tab, AIII.3.counts.tab, AIM.1.counts.tab, AIM.2.counts.tab, AIM.3.counts.tab, AIV.1.counts.tab, AIV.2.counts.tab, AIV.3.counts.tab, C.1.counts.tab, C.2.counts.tab, C.3.counts.tab, E.1.counts.tab, E.2.counts.tab, E.3.counts.tab, G.1.counts.tab, H.1.counts.tab, H.2.counts.tab, M1.1.counts.tab, M2.1.counts.tab, M3.1.counts.tab, AV.1.counts.tab, AV.2.counts.tab, AV.4.counts.tab
#Author: Jessica Goodheart, Gigi Minsky and Michael S. Drummond
#Goal: Take the raw fastq files, perform QC checks, filter the adaptor, remove PCR duplicates, trim, and map reads to reference transcriptome. 

import os, re, sys, string

# Use 'find . -maxdepth 3 -name "fastq" > read_files' to identify the folders and read files of interest
file_list = open("read_files", "r")
ref = sys.argv[1]

# Counter to figure out the number of jobs run
count = 0

for line in file_list:
    #pull out the line and strip off the newline and extra stuff
    path = line.strip()
    directory = path.strip(".fastq").split("/")[0]
    name = path.strip(".fastq").split("/")[1]
    spec = path.split(".")[1].split("/")[0]
    print directory
    count += 1
    print "Processing " + name + " now..."

    # GENERATE COMMANDS FOR BASH SCRIPTS
    ######### - FILTER AND TRIM - #########
    # Start quality filter
    qual = "QualFilterFastq.pl -i " + name + ".fastq -m 20 -x 10 -o " + name + "_trimmed.fastq"
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
    align = "gmapper --qv-offset 33 -Q --strata -o 3 -N 1 " + name + "_trimfiltcomplete.fastq " + ref + " > " + name + ".sam"
    # Filter alignments for weak/ambiguous/duplicated alignments
    filt = "SAMFilterByGene.pl -i " + name + ".sam -m 40 -o " + name + "_filtered.sam -p 1"

    os.system("echo \"#!/bin/bash\n\n#SBATCH --job-name=tagseq\n#SBATCH --workdir=/home/goodheart/projects/upper_lip_expression/tag-seq_data/bio_pathway_mike_gigi/" + directory + "/\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=64G\n#SBATCH --time=72:00:00 \n\n" + qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align + "\n" + filt + "\" > " + name + "_job.sh")

    os.system("sbatch " + name + "_job.sh")

file_list.close()
print "Completed processing " + str(count) + " files."

#qual + "\n" + HRF + "\n" + adapter + "\n" + dups + "\n" + tags + "\n" + align + "\n" + filt + "\" > " + name + "_job.sh")


