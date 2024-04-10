#Description: Find BCN one-to-one orthologs and all BCN orthologs (i.e. many-to-many, one-to-many & etc.) across V.tsujii and Skogsbergia sp. transcriptomes. 
#Author: Lisa Yeter Mesrop 
#Goal: Use Orthofinder output to find all BCN one-to-one orthologs and all BCN orthologs for BCN conservation test. 

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

##### import sheets #####

#import the BCN module file 
red_module <- read.csv("050222_red_module_gene_annot_power_8_prefilter_5_3.csv") 

# subset to just include the transcript ids 
red_module_transcript <- subset(red_module, select = "transcript_id")

#change the column to class character for downstream functions 
red_module_transcript <- red_module_transcript %>% mutate(transcript_id = as.character(transcript_id))

#import the Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp file from the Orthofinder output. 
vargulatsujii_v_skogsbergia_ALL_orthologs <- read.csv("Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp.csv", header = TRUE)

##### use the function to find the number of BCN orthologs #####

# function to check for matches and subset the vargulatsujii_v_skogsbergia_ALL_orthologs df by the matched BCN transcripts. NOTE: when using the check_and subset function, make sure to change the column names accordingly 

check_and_subset <- function(df1, df2) {
  matched_rows_list <- list()
  for (i in 1:nrow(df1)) {
    char_row <- df1[i, "transcript_id"]
    matched_rows <- df2[str_detect(df2$Vargula_tsujii_cdhit_95.fasta.transdecoder, paste0(char_row, collapse = "|")), , drop = FALSE]
    if (nrow(matched_rows) > 0) {
      matched_rows_list[[i]] <- matched_rows
    }
  }
  if (length(matched_rows_list) > 0) {
    matched_df <- bind_rows(matched_rows_list)
    return(matched_df)
  } else {
    return(NULL)
  }
}

# use the function to find all BCN orthologs
BCN_tsujii_and_skogs_orthogroups <- check_and_subset(red_module_transcript, vargulatsujii_v_skogsbergia_ALL_orthologs)

#take a peek 
head(BCN_tsujii_and_skogs_orthogroups)

##### find all BCN orthologs #####

# take the subsetted df and change the Skogsbergia column to a character 
BCN_tsujii_and_skogs_orthogroups <- BCN_tsujii_and_skogs_orthogroups %>% mutate(Skogsbergia_sp = as.character(Skogsbergia_sp))

# create a character vector to hold Skogsbergia sp transcripts
BCN_ALL_orthologs_skogs_ids <- BCN_tsujii_and_skogs_orthogroups$Skogsbergia_sp

#unlist to separate the transcript ids 
BCN_ALL_orthologs_skogs_ids_unlist  <- unlist(strsplit(BCN_ALL_orthologs_skogs_ids,","))

#remove any trailing empty spaces 
BCN_ALL_orthologs_skogs_ids_unlist <-trimws(BCN_ALL_orthologs_skogs_ids_unlist)

#remove .p* at the end of the transcripts 
BCN_ALL_orthologs_skogs_ids_unlist_removep <- BCN_ALL_orthologs_skogs_ids_unlist %>% str_replace("(.p1)", "") %>% str_replace("(.p2)", "")

#make unique 
BCN_ALL_orthologs_skogs_ids_unlist_removep_unique <- unique(BCN_ALL_orthologs_skogs_ids_unlist_removep)

#save 
write.csv(BCN_ALL_orthologs_skogs_ids_unlist_removep_unique, file = "BCN_ALL_orthologs_skogs_ids_unlist_removep_unique.csv")

##### find all BCN one-to-one orthologs #####

# subset the BCN_tsujii_and_skogs_orthogroups to one-to-orthologs. To perform this, remove any rows in V.tsujii column and Skogsbergia sp. column that have multiple orthologs. 

BCN_one_to_one_orthologs_skogs <- BCN_tsujii_and_skogs_orthogroups %>% filter(!grepl(',', Vargula_tsujii_cdhit_95.fasta.transdecoder))
BCN_one_to_one_orthologs_skogs <- BCN_one_to_one_orthologs_skogs %>% filter(!grepl(',', Skogsbergia_sp))

# create a character vector to hold Skogsbergia sp. transcripts
BCN_one_to_one_orthologs_skogs_ids <- BCN_one_to_one_orthologs_skogs$Skogsbergia_sp

# unlist to separate the transcript ids 
BCN_one_to_one_orthologs_skogs_ids_unlist  <- unlist(strsplit(   BCN_one_to_one_orthologs_skogs_ids,","))

# remove any trailing empty spaces 
BCN_one_to_one_orthologs_skogs_ids_unlist <-trimws(BCN_one_to_one_orthologs_skogs_ids_unlist)

# remove .p* at the end of the transcripts 
BCN_one_to_one_orthologs_skogs_ids_unlist_removep <- BCN_one_to_one_orthologs_skogs_ids_unlist %>% str_replace("(.p1)", "") %>% str_replace("(.p2)", "") 

#make unique 
BCN_one_to_one_orthologs_skogs_ids_unlist_removep_unique <- unique(BCN_one_to_one_orthologs_skogs_ids_unlist_removep)

#save
write.csv(BCN_one_to_one_orthologs_skogs_ids_unlist_removep_unique, file = "BCN_one_to_one_orthologs_skogs_ids_unlist_removep_unique.csv")









