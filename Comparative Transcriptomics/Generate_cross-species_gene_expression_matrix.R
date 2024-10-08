#Description: Generate cross-species gene expression matrix for PCA. 
#Author: Lisa Yeter Mesrop
#Goal: Use the individual count expression matrices for V.tsujii and Skogsbergia sp. datasets (upper lip, compound eye, gut) and output from Orthofinder (i.e. one-to-one orthologs) to generate cross-species expression matrix. 

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

#first start with V.tsujii 
vargulatsujii_v_skogsbergia_ALL_orthologs_v1 <- vargulatsujii_v_skogsbergia_ALL_orthologs %>% filter(!grepl(',', Vargula_tsujii_cdhit_95.fasta.transdecoder))

#do the same for Skogsbergia sp. 
vargulatsujii_v_skogsbergia_ALL_orthologs_v2 <- vargulatsujii_v_skogsbergia_ALL_orthologs_v1 %>% filter(!grepl(',', Skogsbergia_sp))

#take a quick peek at the df which now has expression counts of one-to-one orthologs across V.tsujii and Skogsbergia sp. transcriptomes 
head(All_one_to_one_orthlogs_vtsujii_skogs)

##### clean up the df All_one_to_one_orthlogs_vtsujii_skogs #####

#remove extra column and clean up 
All_one_to_one_orthlogs_vtsujii_skogs$X <- NULL

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

##### now identify all one-to-one orthologs present in each expression matrix #####

#start with V.tsujii expression matrix 
vtsujii_one_to_one_ortholog <- subset(vtsujii_up_eye_gut_counts, X %in% All_one_to_one_orthlogs_vtsujii_skogs$tsujii)

#change column name to tsujii 
colnames(vtsujii_one_to_one_ortholog)[1] <- "tsujii"

#next with Skogsbergia sp. expression matrix 
skogs_one_to_one_ortholog <- subset(skogs_up_eye_gut_counts, X %in% All_one_to_one_orthlogs_vtsujii_skogs$skogs)

#change column name to skogs
colnames(skogs_one_to_one_ortholog)[16] <- "skogs"

#add the orthogroup number to each subsetted matrix
skogs_filtered_join  <- plyr::join(skogs_one_to_one_ortholog , All_one_to_one_orthlogs_vtsujii_skogs, by = "skogs", type ="left", match = "all")
vtsujii_filtered_join  <- plyr::join(vtsujii_one_to_one_ortholog , All_one_to_one_orthlogs_vtsujii_skogs, by = "tsujii", type ="left", match = "all")

#take the largest expression matrix and subset it by the smaller expression matrix 
vtsujii_skogs_joined_one_to_one  <- plyr::join(vtsujii_filtered_join ,skogs_filtered_join, by = "tsujii", type ="left", match = "all")

#omit any rows that don't have expression counts in both V.tsujii and Skogsbergia sp. expression matrices
vtsujii_skogs_joined_one_to_one_na_omit <- na.omit(vtsujii_skogs_joined_one_to_one)

#remove any duplicated columns from the plyr::join function 
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders <- vtsujii_skogs_joined_one_to_one_na_omit %>% dplyr::select(-c(17, 18))

#In some cases, one orthogroup can have multiple one-to-one orthologs if the gene duplication occurred before the divergence of the two species. To deal with this, we need to generate unique names for those orthogroup numbers in column Orthogroup. 
#Note: An alternative approach is to combine the expression of multiple one-to-one orthologs within an orthogroup to calculate an average count per orthogroup for each species.
#convert the Orthogroup column from factor to character 
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders <- vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders %>%
  mutate(Orthogroup = as.character(Orthogroup))

#identify duplicated rows based on the column Orthogroup
duplicated_rows <- duplicated(vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders$Orthogroup)

#function to add consecutive numbers to duplicate orthogroup numbers in column Orthogroup to make them unique. 
add_consecutive_num_to_dup_ortho <- function(x) {
  if (any(duplicated(x))) {
    seq_num <- seq_along(x)
    x[duplicated(x)] <- paste0(x[duplicated(x)], "_duplicated", seq_num[duplicated(x)])
  }
  return(x)
}

#use the function to add consecutive numbers to duplicated orthogroup numbers in the column Orthogroup
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders <-vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders %>%
  mutate(Orthogroup = add_consecutive_num_to_dup_ortho(Orthogroup))

#make the column Orthogroup the rownames
row.names(vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders) <- vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders$Orthogroup

#clean up df 
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders$tsujii <- NULL
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders$skogs <- NULL
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders$Orthogroup <- NULL

#final df 
vtsujii_skogs_joined_one_to_one_na_omit_cleanheaders

















