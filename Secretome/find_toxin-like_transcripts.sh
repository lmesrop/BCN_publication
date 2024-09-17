#Description: Identify putative toxin-like genes in the putative secretomes of tissue types and BCN. 
#Author: Lisa Yeter Mesrop 
#Goal: Genes were identified as putative toxin-like if they contained a domain associated with known toxin protein families.

#load libraries 
library(matrixStats)
library(dplyr)
library(data.table)
library(ggraph) 
library(graphlayouts)
library(igraph)
library(readxl)
library(plyr)
library(tidyverse)

#read in the putative secretomes for significantly upregulated genes of upper lips, compound eyes, and guts of Vargula tsujii and Skogsbergia sp. Read in the putative secretome for the BCN. All dataframes include annotations. 

secretome_Skogs_sigfig_unique_UpperLip
secretome_Skogs_sigfig_unique_comEye
secretome_Skogs_sigfig_unique_Gut
secretome_Vtsujii_sigfig_unique_BioUpperLip
secretome_Vtsujii_sigfig_unique_comEye
secretome_Vtsujii_sigfig_unique_Gut
BCN_annotated_SignalPep_NoTransMem

Pfam_toxin_domains <- data.frame(pfam_ids = c("PF16470", "PF13688", "PF06369", "PF01549", "PF01401",
                                              "PF00450", "PF00431", "PF00188", "PF00095", "PF00089", "PF00061","PF08212", "PF00059", "PF00014" , "PF05826", "PF16845", 
                                              "PF00031", "PF00704", "serpin", "Serpin", "Sphingomyelinase C", "sphingomyelinase C", "phospholipase C", "Phospholipase C", "sphingomyelin phosphodiesterase", "Omega-agatoxin", "omega-agatoxin"))


# function to find matches in the dataset of interest
toxin_matches <- function(toxin_ids, dataset) {
  matches <- grep(toxin_ids, dataset, ignore.case = TRUE)
  if (length(matches) > 0) {
    return(matches)
  } else {
    return(NA)
  }
}

#### Vargula tsujii ####

### bioluminescent upper lip  ### 

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Vtsujii_sigfig_unique_BioUpperLip$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
subset_Vtsujii_Bio_Upper_Lip_pfam <- secretome_Vtsujii_sigfig_unique_BioUpperLip[toxin_matched_indices, ]

#rm duplicates 

subset_Vtsujii_Bio_Upper_Lip_pfam_df <- subset_Vtsujii_Bio_Upper_Lip_pfam %>% distinct(transcript_id)
nrow(subset_Vtsujii_Bio_Upper_Lip_pfam_df)

### compound eye ###

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Vtsujii_sigfig_unique_comEye$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
subset_Vtsujii_eye_pfam <- secretome_Vtsujii_sigfig_unique_comEye[toxin_matched_indices, ]

#rm duplicates 

subset_Vtsujii_eye_pfam_df <- subset_Vtsujii_eye_pfam %>% distinct(transcript_id)
nrow(subset_Vtsujii_eye_pfam_df) 

### gut ###

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Vtsujii_sigfig_unique_Gut$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
subset_Vtsujii_gut_pfam <-secretome_Vtsujii_sigfig_unique_Gut[toxin_matched_indices, ]

#rm duplicates 
subset_Vtsujii_gut_pfam_df <- subset_Vtsujii_gut_pfam %>% distinct(transcript_id)
nrow(subset_Vtsujii_gut_pfam_df) 

### BCN ###

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, BCN_annotated_SignalPep_NoTransMem$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
subset_Vtsujii_BCN_pfam <-BCN_annotated_SignalPep_NoTransMem[toxin_matched_indices, ]

#rm duplicates 

subset_Vtsujii_BCN_pfam_df <- subset_Vtsujii_BCN_pfam %>% distinct(transcript_id)
nrow(subset_Vtsujii_BCN_pfam_df) 

#determine how many of subset_Vtsujii_BCN_pfam_df have high connectivity (MM>0.8) in the BCN

subset_Vtsujii_BCN_pfam_df$transcript_id[subset_Vtsujii_BCN_pfam_df$transcript_id %in% top_datKME_ADJ1_.8_desc_BCN$transcript_id]


#### Skogsbergia tissue ####

### upper lip ###

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Skogs_sigfig_unique_UpperLip$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
subset_secretome_Skogs_sigfig_unique_UpperLip_pfam <- secretome_Skogs_sigfig_unique_UpperLip[toxin_matched_indices, ]

#rm duplicates 

subset_secretome_Skogs_sigfig_unique_UpperLip_pfam_df <- subset_secretome_Skogs_sigfig_unique_UpperLip_pfam %>% distinct(transcript_id)
nrow(subset_secretome_Skogs_sigfig_unique_UpperLip_pfam_df) 

### compound eye ###

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Skogs_sigfig_unique_comEye$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
secretome_Skogs_sigfig_unique_comEye_pfam <- secretome_Skogs_sigfig_unique_comEye[toxin_matched_indices, ]

#rm duplicates 

secretome_Skogs_sigfig_unique_comEye_pfam_df <- secretome_Skogs_sigfig_unique_comEye_pfam %>% distinct(transcript_id)
nrow(secretome_Skogs_sigfig_unique_comEye_pfam_df) 


### gut ### 

# apply function to find matches and subset dataframe
toxin_matched_list <- lapply(Pfam_toxin_domains$pfam_ids, function(x) toxin_matches(x, secretome_Skogs_sigfig_unique_Gut$Pfam))
toxin_matched_indices <- unlist(toxin_matched_list[!is.na(toxin_matched_list)])

# subset df2 based on matched indices
secretome_Skogs_sigfig_unique_Gut_pfam <-secretome_Skogs_sigfig_unique_Gut[toxin_matched_indices, ]

# rm duplicates 

secretome_Skogs_sigfig_unique_Gut_pfam_df <- secretome_Skogs_sigfig_unique_Gut_pfam %>% distinct(transcript_id)
nrow(secretome_Skogs_sigfig_unique_Gut_pfam_df) 






