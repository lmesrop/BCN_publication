#!/bin/sh

#SBATCH --job-name=Skogs_trinity
#SBATCH --workdir=/home/lmesrop/upper_lip_expression_project
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=108:00:00

source /home/lmesrop/.bash_profile

Trinity --seqType fq --max_memory 96G --left Sk_1_1.fq.gz --right Sk_1_2.fq.gz --CPU 16 --trimmomatic

