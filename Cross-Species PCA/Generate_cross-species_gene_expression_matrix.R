#Description: Generate cross-species gene expression matrix for downstream analyses. 
#Author: Lisa Yeter Mesrop
#Goal: 

#load libraries 
library(tidyverse) 
library(edgeR) 
library(matrixStats)
library(DESeq2)
library(plyr)
library(dplyr)
library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(repr)
library(topGO)
library(reshape2)
library(scales)

##### Import all orthologs btw V.tsujii and Skogsbergia sp. from the Orthofinder output #####

vargulatsujii_v_skogsbergia_ALL_orthologs <- read.csv("Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp.csv", header = TRUE)

#take a quick peek 
head(vargulatsujii_v_skogsbergia_ALL_orthologs)

##### Subset to include only one-to-one orthologs #####

#find all orthogroups that have one-to-one orthologs across the two columns (i.e. first column V.tsujii orthologs and second column Skogsbergia sp.)  
#remove any rows that have multiple orthologs in each column 

#first start with V.tsujii 
vargulatsujii_v_skogsbergia_ALL_orthologs_v1 <- vargulatsujii_v_skogsbergia_ALL_orthologs %>% filter(!grepl(',', Vargula_tsujii_cdhit_95.fasta.transdecoder))

#do the same for Skogsbergia sp. 
vargulatsujii_v_skogsbergia_ALL_orthologs_v2 <- vargulatsujii_v_skogsbergia_ALL_orthologs_v1 %>% filter(!grepl(',', Skogsbergia_sp))

#take a quick peek at the df which now has expression counts of one-to-one orthologs across V.tsujii and Skogsbergia sp. transcriptomes 
head(All_one_to_one_orthlogs_vtsujii_skogs)

##### clean up the df All_one_to_one_orthlogs_vtsujii_skogs #####

#remove extra column
All_one_to_one_orthlogs_vtsujii_skogs$X <- NULL

#remove .p1 and .p2 at the end of each transcript ID. TIP: To avoid multiple open reading frames for one transcript from TransDecoder output, use the flag --single_best_orf. 

All_one_to_one_orthlogs_vtsujii_skogs_v1 <- All_one_to_one_orthlogs_vtsujii_skogs %>% separate(Skogsbergia_sp, c("skogs", "extra"), sep ="\\.p") 
All_one_to_one_orthlogs_vtsujii_skogs_v1$extra <- NULL
All_one_to_one_orthlogs_vtsujii_skogs_v2 <- All_one_to_one_orthlogs_vtsujii_skogs_v1 %>% separate(Vargula_tsujii_cdhit_95.fasta.transdecoder, c("tsujii", "extra"), sep ="\\.p") 
All_one_to_one_orthlogs_vtsujii_skogs_v2$extra <- NULL
All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup <- All_one_to_one_orthlogs_vtsujii_skogs_v2 %>% distinct(skogs, .keep_all = TRUE)
All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups <- All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup %>% distinct(tsujii, .keep_all = TRUE)  

#cleaned up df 
All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups

##### read in count matrices for V.tsujii and Skogsbergia sp. datasets (upper lip, compound eye and gut) #####

#read in the Skogsbergia sp. dataset 
skogs_up_eye_gut_counts <- read.delim("skogs_fasta90_isoform_combined.tab", header = TRUE, sep = "\t", quote = "")

#clean up the df 
skogs_up_eye_gut_counts$X.1 <- NULL 

#Keep the X column for downstream analyses 

#read in the V.tsujii dataset 
vtsujii_up_eye_gut_counts <- read.delim("combined_vargula_tsujii.tab", header = TRUE, sep = "\t", quote = "")

#clean up the df 
vtsujii_up_eye_gut_counts$X.1 <- NULL

#make X into rownames 
row.names(vtsujii_up_eye_gut_counts) <- vtsujii_up_eye_gut_counts$X

#keep the X column for downstream analyses 


##### now find all the one-to-one orthologs in each expression matrix #####

#start with V.tsujii 
vtsujii_one_to_one_ortholog <- subset(vtsujii_up_eye_gut_counts, X %in% All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups$tsujii)

#change column names 
colnames(vtsujii_one_to_one_ortholog)[1] <- "tsujii"

#next with Skogsbergia sp. 
skogs_one_to_one_ortholog <- subset(skogs_up_eye_gut_counts, X %in% All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups$skogs)

#change column names 
colnames(skogs_one_to_one_ortholog)[16] <- "skogs"

#add the orthogroup information to each subsetted matrix
skogs_filtered_join  <- plyr::join(skogs_one_to_one_ortholog , All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups, by = "skogs", type ="left", match = "all")
vtsujii_filtered_join  <- plyr::join(vtsujii_one_to_one_ortholog , All_one_to_one_orthlogs_vtsujii_skogs_v2_rmskogsdup_rmtsujii_dups, by = "tsujii", type ="left", match = "all")

#take the largest expression matrix and subset it by the smaller matrix 
vtsujii_skogs_joined_one_to_one  <- plyr::join(vtsujii_filtered_join ,skogs_filtered_join, by = "tsujii", type ="left", match = "all")

#omit any rows that don't have expression counts in both V.tsujii and Skogsbergia sp. 
vtsujii_skogs_joined_one_to_one_na_omit <- na.omit(vtsujii_skogs_joined_one_to_one)

#the plyr::join function creates duplicated columns. Remove the duplicated columns. 
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders <- vtsujii_skogs_joined_one_to_one_na_omit %>% dplyr::select(-c(17, 18))

#remove duplicated orthogroup numbers ?
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders_non_dup <-  vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders %>% distinct(Orthogroup, .keep_all = TRUE)

#





















