#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --time=208:00:00
#SBATCH --mail-user=lmesrop@ucsb.edu
#SBATCH --job-name=transdecoder

cd $SLURM_SUBMIT_DIR

TransDecoder.LongOrfs -t /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform.fasta 
TransDecoder.Predict -t /home/lmesrop/skogs_belize_transcriptome/skogs_cdhit_90_longest_isoform_fasta_transdecoder_bowtie2_050922/skogs_cdhit_90_longestisoform.fasta

