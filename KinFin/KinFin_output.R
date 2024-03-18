#Title: Downstream analyses of KinFin output for Vargula tsujii 
#Author: Lisa Yeter Mesrop 
#Description: Determine the distribution of Vargula tsujii-specific genes across different taxonomic origins. 

#load libraries 
library(tidyverse) 
library(edgeR)
library(matrixStats)
library(DESeq2)
library(data.table)
library(ggraph) 
library(graphlayouts)
library(igraph)
library(topGO)
library(readxl)
library(dplyr)
library(plyr)
library(ggvenn) 

##### upload KinFin outputs and dataset dataframes of interest #####

# upload KinFin outputs for V.tsujii (without nested proteins)

Vargula_tsujii_cdhit_95_proteins
arthro_only_proteins
lum_only_proteins
lux_only_proteins
ostra_only_proteins

#remove .p1 or .p2 
Vargula_tsujii_cdhit_95_proteins_removep <- Vargula_tsujii_cdhit_95_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
ostra_only_proteins_removep <- ostra_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
arthro_only_proteins_removep <-  arthro_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lum_only_proteins_removep <- lum_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lux_only_proteins_removep <- lux_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))

#dataframes for significantly upregulated genes of the luminous upper lip, compound eye and gut
#luminous upper lip 
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene
#compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene
#gut 
sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene

#dataframes for co-expressed BCN genes 
#BCN
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id
#gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id


##### determine the number of ORFs in each dataset of interest #####

#create a dataframe for the total number of transcripts that have predicted ORFs (and used for OrthoFinder and KinFin analyses). 
#combined arthro_proteins.txt and Vargula_tsujii_cdhit_95_proteins.txt from KinFin output. 
#remove duplicated transcripts. 

ALL_arthro_and_vtsujii_specific_proteins_rm <- ALL_arthro_and_vtsujii_specific_proteins %>% distinct()
nrow(ALL_arthro_and_vtsujii_specific_proteins_rm) 
ALL_arthro_and_vtsujii_specific_proteins_rm_removep <- ALL_arthro_and_vtsujii_specific_proteins_rm %>% separate(...1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))

#for each dataset, determine the number of ORFS for each taxonomic level of interest. 
#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id]

#bio upper lip 
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id]

#compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id]

#gut
gut_orfs_total <- sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id]

#gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id]

#all expressed genes in the DGE dataset 
#DGE dataset includes (bio upper lip, compound eye and gut) with a min of 2 counts or more in greater than 5 biological replicates 
All_expressed_transcripts_arthro_vtsujii_specific_prefiltered <- dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids %in% ALL_arthro_and_vtsujii_specific_proteins_rm_removep$transcript_id] 

##### determine the number of Vargula tsujii specific genes across each taxonomic origin (from KinFin output) #####

### V.tsujii ###

#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#Bio upper lip 
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#Compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#Gut
sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#Gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates)
All_expressed_transcripts_vtsujii_terminal_orphan_min2counts <- dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]


### Ostracoda ###

#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% ostra_only_proteins_removep$transcript_id ]

#bio upper lip 
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% ostra_only_proteins_removep$transcript_id ]

#compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% ostra_only_proteins_removep$transcript_id ]

#gut
sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% ostra_only_proteins_removep$transcript_id ]

#gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% ostra_only_proteins_removep$transcript_id ]

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates)
All_expressed_trascripts_ostrac_strict_min5 <- dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids  %in% ostra_only_proteins_removep$transcript_id]

### Arthropoda ###

#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% arthro_only_proteins_removep$transcript_id]

#bio upper lip 
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% arthro_only_proteins_removep$transcript_id]

#compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% arthro_only_proteins_removep$transcript_id]

#gut
arthropod_sig_gut <- sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% arthro_only_proteins_removep$transcript_id]

#gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% arthro_only_proteins_removep$transcript_id]

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates)
ALL_expressed_transcripts_arthro_min5 <- dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids  %in% arthro_only_proteins_removep$transcript_id]

### Luminini ###

#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% lum_only_proteins_removep$transcript_id]

#bio upper lip 
LUM_DE_BIOUPPER <- sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% lum_only_proteins_removep$transcript_id]

#compound eye
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% lum_only_proteins_removep$transcript_id]

#gut
sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% lum_only_proteins_removep$transcript_id]

#gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% lum_only_proteins_removep$transcript_id]


#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates or more)
dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids %in% lum_only_proteins_removep$transcript_id]
                                                                 

### Luxorina ###

#BCN 
X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id[X050222_red_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% lux_only_proteins_removep$transcript_id]

#bio upper lip
sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% lux_only_proteins_removep$transcript_id]

#compound eye 
sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene[sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene %in% lux_only_proteins_removep$transcript_id]

#gut
sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% lux_only_proteins_removep$transcript_id]

#Gut module 
vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id[vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3$transcript_id %in% lux_only_proteins_removep$transcript_id]

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates or more)
dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids[dds_count_table_organ_level_prefiltered_strict_assay_ids$transript_ids %in% lux_only_proteins_removep$transcript_id]


##### determine the number of Vargula tsujii specific genes across each taxonomic origin (from KinFin output) UNIQUELY expressed in each dataset #####

#create dfs with transcript ids for each of the significantly upregulated genes of bio upper lip, compound eye and gut

lightorgan_de <- sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene

eye_de <- sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene

gut_de <- sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene

#generate a venn diagram to visualize shared significantly upregulated genes across tissue types. 
unique_venn_list <- list(
  Bio_Upper_Lip =  lightorgan_de , 
  Gut = gut_de,
  Compound_Eye = eye_de
)

ggvenn(
  unique_venn_list, 
  fill_color = c("#3DB3C2", "#71C23D", "#8E3DC2"),
  stroke_size = .7, set_name_size = 4.5
)


#find the genes uniquely upregulated in the bio upper lip
lightorgan_unique_de  <- sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene[!sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene %in% sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene]
lightorgan_unique_de_df <- as.data.frame(lightorgan_unique_de)
nrow(lightorgan_unique_de_df)
#476

#find the genes uniquely upregulated in the gut 
gut_unique_de_v1   <- sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene[!sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene %in% sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122$gene]
#2216 
#now subtract the compound eye from gut_unique_de_v1
gut_unique_de <- gut_unique_de_v1[!gut_unique_de_v1 %in% sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122$gene]
gut_unique_de_df <- as.data.frame(gut_unique_de)
nrow(gut_unique_de_df)
#2206

#find the genes uniquely upregulated in the compound eye 
eye_unique_de <- eye_de[!eye_de %in% sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122$gene]
eye_unique_de_df <- as.data.frame(eye_unique_de)
nrow(eye_unique_de_df)
#350 

#determine the number of significantly upregulated genes uniquely expressed in the bio upper lip 

#vargula tsujii 
lightorgan_unique_de_df$lightorgan_unique_de[lightorgan_unique_de_df$lightorgan_unique_de %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#ostracoda 
lightorgan_unique_de_df$lightorgan_unique_de[lightorgan_unique_de_df$lightorgan_unique_de %in% ostra_only_proteins_removep$transcript_id]

#arthropoda 
lightorgan_unique_de_df$lightorgan_unique_de[lightorgan_unique_de_df$lightorgan_unique_de %in% arthro_only_proteins_removep$transcript_id]

#luminini 
lightorgan_unique_de_df$lightorgan_unique_de[lightorgan_unique_de_df$lightorgan_unique_de %in% lum_only_proteins_removep$transcript_id]

#luxorina
lightorgan_unique_de_df$lightorgan_unique_de[lightorgan_unique_de_df$lightorgan_unique_de %in% lux_only_proteins_removep$transcript_id]


#determine the number of significantly upregulated genes uniquely expressed in the compound eye 

#vargula tsujii 
eye_unique_de_df$eye_unique_de[eye_unique_de_df$eye_unique_de %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#ostracoda 
eye_unique_de_df$eye_unique_de[eye_unique_de_df$eye_unique_de %in% ostra_only_proteins_removep$transcript_id]

#arthropoda 
eye_unique_de_df$eye_unique_de[eye_unique_de_df$eye_unique_de %in% arthro_only_proteins_removep$transcript_id]

#luminini 
eye_unique_de_df$eye_unique_de[eye_unique_de_df$eye_unique_de %in% lum_only_proteins_removep$transcript_id]

#luxorina 
eye_unique_de_df$eye_unique_de[eye_unique_de_df$eye_unique_de %in% lux_only_proteins_removep$transcript_id]


#determine the number of significantly upregulated genes uniquely expressed in the gut 

#vargula tsujii 
gut_unique_de_df$gut_unique_de[gut_unique_de_df$gut_unique_de %in% Vargula_tsujii_cdhit_95_proteins_removep$transcript_id]

#ostracoda 
gut_unique_de_df$gut_unique_de[gut_unique_de_df$gut_unique_de %in% ostra_only_proteins_removep$transcript_id]

#arthropoda 
gut_unique_de_df$gut_unique_de[gut_unique_de_df$gut_unique_de %in% arthro_only_proteins_removep$transcript_id]

#luminini 
gut_unique_de_df$gut_unique_de[gut_unique_de_df$gut_unique_de %in% lum_only_proteins_removep$transcript_id]

#luxorina 
gut_unique_de_df$gut_unique_de[gut_unique_de_df$gut_unique_de %in% lux_only_proteins_removep$transcript_id]



