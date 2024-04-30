
# Description: Generate GO-Figure! semantic plots 
# Package: GO-Figure! 
# Goal: Generate GO semantic plots using outputs from TopGO enrichment analyses performed in Jupyter Notebook.  

### Script 1

Go-Figure! analysis for both species, Vargula tsujii, and Skogsbergia sp. using DE genes (significantly upregulated) for both luminous and non-luminous upper lips. 

python /home/lmesrop/programs/GO-Figure/gofigure.py -j standard-plus -i GO_figure_bothspecies_df_subset.tsv -o gofigure_directory_DE_bothspecies_si0.5 -si 0.5 -c user -m 50 -e 100 -q pdf

### Script 2 

Go-Figure! analysis for Skogsbergia sp. using DE genes (significantly upregulated) for non-luminous upper lip. 

python /home/lmesrop/programs/GO-Figure/gofigure.py -j topgo -i GO_figure_DE_Skogsbergia_sp_df.tsv -o gofigure_directory_DE_skogsbergia_sp_si0.5 -si 0.5 -m 25 -e 100 -q pdf


### Script 3

Go-Figure! analysis for Vargula tsujii using DE genes (significantly upregulated) for luminous upper lip. 

python /home/lmesrop/programs/GO-Figure/gofigure.py -j topgo -i GO_figure_DE_bio_upper_lip_Vargulatsujii_df_nodes61.tsv -o GO_BioUL_gofigure_directory_si0.5 -e 100 -m 100 -si 0.5 -q pdf



### Script 4 

GO-Figure! analysis for all co-expressed genes of the BCN (Vargula tsujii). 

python /home/lmesrop/programs/GO-Figure/gofigure.py -j topgo -i GO_figure_BCN_df.tsv -o gofigure_directory_BCN_si0.1_color -si 0.1 -m 50 -e 100 -q pdf -p Purples



