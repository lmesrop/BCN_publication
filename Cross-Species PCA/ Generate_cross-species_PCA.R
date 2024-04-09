#Description: Generate the cross-species PCA plot for Figure 3 
#Author: Lisa Yeter Mesrop 
#Goal: Take the cross-species expression matrix (i.e. one-to-one orthologs) and filter, log2 transform and remove batch effect caused by multiple species. 

#load libraries 

library(tidyverse) 
library(DESeq2)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(plyr)
library(scales)
library(sva)
library(RColorBrewer)
library(pheatmap)
library(PCAtools)
library(edgeR)
library(DESeq2)

##### import dataframes  #####

## import the cross-species expression matrix 
vtsujii_skogs_joined_one_to_one_na_omit_df_non_dup

## make the metasheet df
meta <- data.frame(row.names = colnames(vtsujii_skogs_joined_one_to_one_na_omit_df_non_dup)) 
sample_name = c("Tsujii_Upper_lip", "Tsujii_Eye", "Tsujii_Gut", "Tsujii_Upper_lip", "Tsujii_Eye", "Tsujii_Gut", "Tsujii_Upper_lip", "Tsujii_Eye", "Tsujii_Gut", "Tsujii_Upper_lip", "Tsujii_Eye", "Tsujii_Gut","Tsujii_Upper_lip", "Tsujii_Eye", "Tsujii_Gut","Skogs_Upper_lip", "Skogs_Eye", "Skogs_Gut", "Skogs_Upper_lip", "Skogs_Eye", "Skogs_Gut", "Skogs_Upper_lip", "Skogs_Eye", "Skogs_Gut", "Skogs_Upper_lip", "Skogs_Eye", "Skogs_Gut","Skogs_Upper_lip", "Skogs_Eye", "Skogs_Gut" )
meta$sample_name <- sample_name 
meta$names <- rownames(meta)
rownames(meta) <- NULL
species = c("Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii","Vargula_tsujii", "Vargula_tsujii", "Vargula_tsujii",
                "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp","Skogsbergia_sp", "Skogsbergia_sp", "Skogsbergia_sp" )
meta$species <- species
tissue <- c("upper_lip", "compound_eye", "gut", "upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut",
          "upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut","upper_lip", "compound_eye", "gut")
meta$tissue <- tissue 

##### prep cross-species expression matrix for PCA #####

## use DESeq2 to remove expression counts with less than 5 counts in more than 3 samples 
dds_count_table <- DESeqDataSetFromMatrix(countData = vtsujii_skogs_joined_one_to_one_na_omit_df_non_dup, colData = meta, design = ~sample_name)
dds_merged_table_prefiltered <- dds_count_table[rowSums(counts(dds_count_table) >= 5) >=3,];
nrow(dds_merged_table_prefiltered)
dds_merged_table_prefiltered_assay <- assay(dds_merged_table_prefiltered)

## use DESeq2 to perform log transformation 
 
vtsujii_skogs_cross_sp_exp_log2 = log2(dds_merged_table_prefiltered_assay + 1e-5) 
row.names(vtsujii_skogs_cross_sp_exp_log2) = row.names(dds_merged_table_prefiltered_assay)
colnames(vtsujii_skogs_cross_sp_exp_log2) = colnames(dds_merged_table_prefiltered_assay)

## remove batch effect with ComBat 

batch = meta$species
vtsujii_skogs_combat = model.matrix(~1, data=meta)
vtsujii_skogs_combat_batch_rm = ComBat(dat=vtsujii_skogs_cross_sp_exp_log2, batch=batch, mod=vtsujii_skogs_combat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)
vtsujii_skogs_combat_batch_rm_Summ_Exp <- SummarizedExperiment(vtsujii_skogs_combat_batch_rm - rowMeans(vtsujii_skogs_combat_batch_rm),colData=meta)
pca_vtsujii_vtsujii_skogs_combat_batch_rm_Summ_Exp <- plotPCA(DESeqTransform(vtsujii_skogs_combat_batch_rm_Summ_Exp), intgroup=c("species", "tissue"), ntop = 50, returnData=TRUE) 
round_percentVar <- round(100 * attr(pca_vtsujii_vtsujii_skogs_combat_batch_rm_Summ_Exp, "percentVar"))
names(pca_vtsujii_vtsujii_skogs_combat_batch_rm_Summ_Exp)[c(4,5)] = c("species", "tissue")

## plot the PCA with ggplot 

ggplot(pca_vtsujii_vtsujii_skogs_combat_batch_rm_Summ_Exp , aes(PC1, PC2, color= tissue, group = tissue ,shape=species)) +
  geom_point(size=3,alpha = 1) +
  scale_shape_manual(values = c(1,17))+
  scale_color_manual(values = c('#C24C3D','#8E3DC2','#E69F00')) +
  xlab(paste0("PC1 (",round_percentVar[1],"%)")) +
  ylab(paste0("PC2 (",round_percentVar[2],"%)")) + 
  coord_fixed() + 
  theme(panel.grid.major = element_line(colour = "gray97",  size = 0.5), panel.grid.minor = element_line(linetype = "dotted"), panel.background = element_rect(fill = NA), 
   legend.key = element_rect(fill = "gray100")) + theme(axis.line = element_line(size = 0.5,linetype = "solid")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  stat_ellipse(type = "norm",geom = "polygon",aes(fill = tissue), alpha = 0.0, level = 0.95,show.legend = F)



