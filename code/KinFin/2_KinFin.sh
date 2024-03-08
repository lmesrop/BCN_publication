
#!/usr/bin/bash

#obtain lists of clade-specific genes

#Adapted from Jessica A. Goodheart (Reference Paper)

#Usage: bash post_kinfin.sh [species] [lux] [lum] [ostra] [arthro] 

#Example: bash post_kinfin.sh Vargula_tsujii_cdhit_95 lux lum ostra arthro

####
# Set up clade names of interest
sp=$1
lx=$2
lm=$3
os=$4
ar=$5

#### Species-Level 
# Get list of species-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' TAXON/TAXON.Vargula_tsujii_cdhit_95.cluster_metrics.txt > "${sp}"_cluster_numbers.txt
awk '{ print $1 }' ../Vtsujii_UnassignedGenes.tsv >> "${sp}"_cluster_numbers.txt

# Get species-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${sp}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${sp}"_clusters.txt
tail -n +2 ../Vtsujii_UnassignedGenes.tsv >> "${sp}"_clusters.txt

# Get lists of species-specific proteins and genes
cat "${sp}"_clusters.txt | tr " " "\n" | grep "NODE_*" | sort | uniq > "${sp}"_proteins.txt
#cat "${sp}"_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${sp}"_genes.txt

#### Luxorina
# Get list of Luxorina-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' luxorina/luxorina."${lx}".cluster_metrics.txt  > "${lx}"_cluster_numbers.txt

# Get Luxorina-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${lx}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${lx}"_clusters.txt

# Get lists of Luxorina-specific proteins and genes
cat "${lx}"_clusters.txt | tr " " "\n" | grep "NODE_*" | sort | uniq > "${lx}"_proteins.txt
comm -23 "${lx}"_proteins.txt "${sp}"_proteins.txt > "${lx}"_only_proteins.txt
#cat "${lx}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${lx}"_only_genes.txt

##### Luminini 
# Get list of Luminini-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' luminini/luminini."${lm}".cluster_metrics.txt > "${lm}"_cluster_numbers.txt

# Get Luminini -specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${lm}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${lm}"_clusters.txt

# Get lists of Luminini -specific proteins and genes
cat "${lm}"_clusters.txt | tr " " "\n" | grep "NODE_*" | sort | uniq > "${lm}"_proteins.txt
comm -23 "${lm}"_proteins.txt "${lx}"_proteins.txt > tmp
comm -23 tmp "${sp}"_proteins.txt > "${lm}"_only_proteins.txt
#cat "${lm}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${lm}"_only_genes.txt

#### Ostracoda 
# Get list of Ostracoda-specific clusters 
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' ostracoda/ostracoda."${os}".cluster_metrics.txt > "${os}"_cluster_numbers.txt

# Get Ostracoda-specific cluster data 
awk 'NR==FNR {a[$1]++; next} $1 in a' "${os}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${os}"_clusters.txt

# Get lists of Ostracoda-specific proteins and genes  
cat "${os}"_clusters.txt | tr " " "\n" | grep "NODE_*" | sort | uniq > "${os}"_proteins.txt
comm -23 "${os}"_proteins.txt "${lm}"_proteins.txt > tmp
comm -23 tmp "${lx}"_proteins.txt > tmp.1
comm -23 tmp.1 "${sp}"_proteins.txt > "${os}"_only_proteins.txt
#cat "${os}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${os}"_only_genes.txt

##### Arthropoda
# Get list of Arthropoda-specific clusters
awk -F"\t" '($2=="present" && $3=="specific") { print $1 }' arthropoda/arthropoda."${ar}".cluster_metrics.txt > "${ar}"_cluster_numbers.txt

# Get Arthropoda-specific cluster data
awk 'NR==FNR {a[$1]++; next} $1 in a' "${ar}"_cluster_numbers.txt ../Orthogroups_edited.txt > "${ar}"_clusters.txt

# Get lists of Arthropoda-specific proteins and genes
cat "${ar}"_clusters.txt | tr " " "\n" | grep "NODE_*" | sort | uniq > "${ar}"_proteins.txt
comm -23 "${ar}"_proteins.txt "${os}"_proteins.txt > tmp
comm -23 tmp "${lm}"_proteins.txt > tmp.1
comm -23 tmp.1 "${lx}"_proteins.txt > tmp.2
comm -23 tmp.2 "${sp}"_proteins.txt > "${ar}"_only_proteins.txt
#cat "${ar}"_only_proteins.txt | sed "s/.t[0-9]*//g" | sort | uniq > "${ar}"_only_genes.txt

rm tmp*

#### Summary statistics 
# Only clade-specific proteins, not including nested proteins
spp=$(cat "${sp}"_proteins.txt | wc -l)
lxp=$(cat "${lx}"_only_proteins.txt | wc -l)
lmp=$(cat "${lm}"_only_proteins.txt | wc -l)
osp=$(cat "${os}"_only_proteins.txt | wc -l)
arp=$(cat "${ar}"_only_proteins.txt | wc -l)

prots=($spp $lxp $lmp $osp $arp)

# ALL clade-specific genes 
sp_all=$(cat "${sp}"_proteins.txt | wc -l)
lx_all=$(cat "${lx}"_proteins.txt | wc -l)
lm_all=$(cat "${lm}"_proteins.txt | wc -l)
os_all=$(cat "${os}"_proteins.txt | wc -l)
ar_all=$(cat "${ar}"_proteins.txt | wc -l)

all_proteins=($sp_all $lx_all $lm_all $os_all $ar_all)

# Create stats file
printf 'clade\tclade_specific_proteins\tall_proteins\n' > kinfin_stats.txt
a=0
for i in "$@"; do
    printf "%s\t%s\t%s\n" "$i" "${prots[a]}" "${all_proteins[a]}" >> kinfin_stats.txt
    ((a+=1))
done
