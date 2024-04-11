#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --time=208:00:00
#SBATCH --mail-user=lmesrop@ucsb.edu
#SBATCH --job-name=cdhit_90 

cd $SLURM_SUBMIT_DIR

cd-hit-est -i /home/lmesrop/skogs_belize_transcriptome/trinity_out_dir/Trinity.fasta -o skogs_cdhit_90 -c 0.90 -n 11 -M 96000 -T 24
