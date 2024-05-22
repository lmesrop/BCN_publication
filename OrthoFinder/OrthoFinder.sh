#!/bin/bash

#SBATCH --nodes=2 --ntasks-per-node=12
#SBATCH --time=990:00:00
#SBATCH --mail-user=lmesrop@ucsb.edu
#SBATCH --job-name=ortho

cd $SLURM_SUBMIT_DIR

orthofinder -f ~/BCN_pub_orthofinder -S diamond -T iqtree -M msa -s ~/BCN_pub_species_orthofinder_iTOL.txt  -y

