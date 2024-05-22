# Description: Steps for KinFin Analysis 
# Author: Adapted from Goodheart et al. (2024) 
# Goal: Prep for KinFin Analysis

#### Step 1: Create config.txt file 

#### Step 2: Copy Orthofinder results into the working directory 

cp Orthogroups.txt /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

cp  SequenceIDs.txt /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

cp SpeciesIDs.txt /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

cp SpeciesTree_rooted.txt /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

cp Orthogroups.tsv /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

cp Orthogroups_UnassignedGenes.tsv /home/lmesrop/wgcna_DE_project_phylogenies/gene_duplication_BCN_phylo_analysis/KinFin_dir/

#### Step 3: Extract V.tsujii-specific genes from Orthogroups_UnassignedGenes.tsv and modify Orthogroups.txt file  

awk -F"\t" '$13 != "" { print $1, $13 }' Orthogroups_UnassignedGenes.tsv > Vtsujii_UnassignedGenes.tsv

sed 's/\t\t/ /g' Orthogroups.tsv | sed 's/\t/ /g' | sed 's/,//g' | sed 's/^M$//g' | sed -r 's/(OG[0-9][0-9][0-9][0-9][0-9][0-9][0-9])\ +/\1\ /g' > Orthogroups_edited.txt

dos2unix Orthogroups_edited.txt

tail -n +2 Orthogroups_edited.txt > Orthogroups_edited.txt.temp

mv Orthogroups_edited.txt.temp Orthogroups_edited.txt

#### Step 4: Run KinFin.sh script 

/home/lmesrop/kinfin/kinfin -g Orthogroups_edited.txt -c KinFin_config_run.txt -s SequenceIDs.txt -p SpeciesIDs.txt -t SpeciesTree_rooted.txt --infer_singletons -o kinfin_results_run


#References: Goodheart, Jessica A et al. “A chromosome-level genome for the nudibranch gastropod Berghia stephanieae helps parse clade-specific gene expression in novel and conserved phenotypes.” BMC biology vol. 22,1 9. 17 Jan. 2024, doi:10.1186/s12915-024-01814-3
