#Description: Downstream analyses of KinFin output for Vargula tsujii 
#Author: Lisa Yeter Mesrop 
#Goal: Determine the distribution of Vargula tsujii-specific genes across different taxonomic origins. 

##### import KinFin outputs and dataset dataframes of interest #####

## import KinFin outputs for V.tsujii (without nested proteins)

Vargula_tsujii_cdhit_95_proteins
arthro_only_proteins
lum_only_proteins
lux_only_proteins
ostra_only_proteins

## clean up dataframes
Vargula_tsujii_cdhit_95_proteins_removep <- Vargula_tsujii_cdhit_95_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
ostra_only_proteins_removep <- ostra_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
arthro_only_proteins_removep <-  arthro_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lum_only_proteins_removep <- lum_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lux_only_proteins_removep <- lux_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))

## import dataframes for significantly upregulated (unique) of the luminous upper lip, compound eye and gut

#luminous upper lip 
df_bio_upper_lip
#compound eye 
df_compound_eye
#gut 
df_gut

#dataframes for co-expressed BCN genes 
#BCN
red_module_gene_annot_power_8_prefilter_5_3$transcript_id

##extract the transcript ids and rename for clarity 
#BCN
BCN_module <- subset(red_module_gene_annot_power_8_prefilter_5_3, select = "transcript_id") 

##create a dataframe for the total number of transcripts that have predicted ORFs (from the OrthoFinder and KinFin analyses). Take the txt files from the KinFin output and combine them into a single dataframe. 
ALL_arthro_and_vtsujii_specific_proteins_rm <- ALL_arthro_and_vtsujii_specific_proteins %>% distinct()
ALL_arthro_and_vtsujii_specific_proteins_rm_removep <- ALL_arthro_and_vtsujii_specific_proteins_rm %>% separate(...1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
#rename for clarity 
ALL_transcripts_with_predictedORFs<- ALL_arthro_and_vtsujii_specific_proteins_rm_removep

#all transcripts with predicted ORFs
ALL_transcripts_with_predictedORFs

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates or more). Imported from Juypter Notebook. 
ALL_DGE_min2counts <- subset(dds_count_table_organ_level_prefiltered_strict_assay_ids, select = "transript_ids")

#### function for finding taxonomically restricted genes for each dataset of interest (e.g. significantly upregulated genes and co-expressed genes)

## first, combine the dataframes from the KinFin output for each taxonomic level (w/o nested proteins)

list1 <- list(arthro_only_proteins_removep, ostra_only_proteins_removep, lum_only_proteins_removep, lux_only_proteins_removep, Vargula_tsujii_cdhit_95_proteins_removep, ALL_transcripts_with_predictedORFs)

names_list_1 <- c("arthro_only_proteins_removep", "ostra_only_proteins_removep", "lum_only_proteins_removep", "lux_only_proteins_removep", "Vargula_tsujii_cdhit_95_proteins_removep", "ALL_transcripts_with_predictedORFs")
  
## rename the lists to the corresponding taxonomic level 
names(list1) <- names_list_1

#function to input the dataset of interest and determine how many are restricted to each taxonomic level of interest and print out the number of transcripts for each taxonomic level 

match_transcripts_to_taxonomic_origin <- function(input_df, match_list) {
  match_chars <- function(row, df) {
    intersect(row, unlist(df))
  }
  matched_dfs <- list()

  for (i in seq_along(match_list)) {
    match_df_name <- names(match_list)[i]  
    match_df <- match_list[[i]]
    
    matched_rows <- purrr::map(input_df, match_chars, match_df)
    
    matched_df <- purrr::map_df(matched_rows, ~ data.frame(matched_characters = .)) 
    
    matched_dfs[[match_df_name]] <- matched_df
  }
  
  return(matched_dfs)
  
  }


#### use the function to find the number of transcripts at each taxonomic level for each dataset ####

### bioluminescent upper lip ###

## call the function 
matched_transcripts_sigfig_up_bio_upper_lip <- match_transcripts_to_taxonomic_origin(df_bio_upper_lip, list1)

## print the output 

for (df_name in names(matched_transcripts_sigfig_up_bio_upper_lip)) {
  df <- matched_transcripts_sigfig_up_bio_upper_lip[[df_name]]
  cat("Dataset: sigfig_up_bio_upper_lip \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_sigfig_up_bio_upper_lip_ARTHRO <-matched_transcripts_sigfig_up_bio_upper_lip[["arthro_only_proteins_removep"]]
matched_sigfig_up_bio_upper_lip_OSTRA <-matched_transcripts_sigfig_up_bio_upper_lip[["ostra_only_proteins_removep"]]
matched_sigfig_up_bio_upper_lip_LUM <- matched_transcripts_sigfig_up_bio_upper_lip[["lum_only_proteins_removep"]]
matched_sigfig_up_bio_upper_lip_LUX <-matched_transcripts_sigfig_up_bio_upper_lip[["lux_only_proteins_removep"]]
matched_sigfig_up_bio_upper_lip_VTSUJII <-matched_transcripts_sigfig_up_bio_upper_lip[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_sigfig_up_bio_upper_lip_ALL_transcripts_with_predictedORFs <-matched_transcripts_sigfig_up_bio_upper_lip[["ALL_transcripts_with_predictedORFs"]]

### compound eye ###

## call the function 
matched_transcripts_sigfig_up_compound_eye <- match_transcripts_to_taxonomic_origin(df_compound_eye, list1)

## print the output 

for (df_name in names(matched_transcripts_sigfig_up_compound_eye)) {
  df <- matched_transcripts_sigfig_up_compound_eye[[df_name]]
  cat("Dataset: sigfig_up_compound_eye \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_sigfig_up_compound_eye_ARTHRO <-matched_transcripts_sigfig_up_compound_eye[["arthro_only_proteins_removep"]]
matched_sigfig_up_compound_eye_OSTRA <-matched_transcripts_sigfig_up_compound_eye[["ostra_only_proteins_removep"]]
matched_sigfig_up_compound_eye_LUM <- matched_transcripts_sigfig_up_compound_eye[["lum_only_proteins_removep"]]
matched_sigfig_up_compound_eye_LUX <-matched_transcripts_sigfig_up_compound_eye[["lux_only_proteins_removep"]]
matched_sigfig_up_compound_eye_VTSUJII <-matched_transcripts_sigfig_up_compound_eye[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_sigfig_up_compound_eye_ALL_transcripts_with_predictedORFs <-matched_transcripts_sigfig_up_compound_eye[["ALL_transcripts_with_predictedORFs"]]


### gut ###

## call the function 
matched_transcripts_sigfig_up_gut <- match_transcripts_to_taxonomic_origin(df_gut, list1)

## print the output 

for (df_name in names(matched_transcripts_sigfig_up_gut)) {
  df <- matched_transcripts_sigfig_up_gut[[df_name]]
  cat("Dataset: sigfig_up_gut \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_sigfig_up_gut_ARTHRO <-matched_transcripts_sigfig_up_gut[["arthro_only_proteins_removep"]]
matched_sigfig_up_gut_OSTRA <-matched_transcripts_sigfig_up_gut[["ostra_only_proteins_removep"]]
matched_sigfig_up_gut_LUM <- matched_transcripts_sigfig_up_gut[["lum_only_proteins_removep"]]
matched_sigfig_up_gut_LUX <-matched_transcripts_sigfig_up_gut[["lux_only_proteins_removep"]]
matched_sigfig_up_gut_VTSUJII <-matched_transcripts_sigfig_up_gut[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_sigfig_up_gut_ALL_transcripts_with_predictedORFs <-matched_transcripts_sigfig_up_gut[["ALL_transcripts_with_predictedORFs"]]

### BCN module dataset ### 

## call the function 
matched_transcripts_BCN_module <- match_transcripts_to_taxonomic_origin(BCN_module, list1)

## print the output 

for (df_name in names(matched_transcripts_BCN_module)) {
  df <- matched_transcripts_BCN_module[[df_name]]
  cat("Dataset: BCN_module \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_BCN_module_ARTHRO <-matched_transcripts_BCN_module[["arthro_only_proteins_removep"]]
matched_BCN_module_OSTRA <-matched_transcripts_BCN_module[["ostra_only_proteins_removep"]]
matched_BCN_module_LUM <- matched_transcripts_BCN_module[["lum_only_proteins_removep"]]
matched_BCN_module_LUX <-matched_transcripts_BCN_module[["lux_only_proteins_removep"]]
matched_BCN_module_VTSUJII <-matched_transcripts_BCN_module[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_BCN_module_ALL_transcripts_with_predictedORFs <-matched_transcripts_BCN_module[["ALL_transcripts_with_predictedORFs"]]

### ALL_DGE_min2counts ###

## call the function 
matched_transcripts_ALL_DGE_min2counts <- match_transcripts_to_taxonomic_origin(ALL_DGE_min2counts, list1)

## print the output 

for (df_name in names(matched_transcripts_ALL_DGE_min2counts)) {
  df <- matched_transcripts_ALL_DGE_min2counts[[df_name]]
  cat("Dataset: ALL_DGE_min2counts \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_ALL_DGE_min2counts_ARTHRO <-matched_transcripts_ALL_DGE_min2counts[["arthro_only_proteins_removep"]]
matched_ALL_DGE_min2counts_OSTRA <-matched_transcripts_ALL_DGE_min2counts[["ostra_only_proteins_removep"]]
matched_ALL_DGE_min2counts_LUM <- matched_transcripts_ALL_DGE_min2counts[["lum_only_proteins_removep"]]
matched_ALL_DGE_min2counts_LUX <-matched_transcripts_ALL_DGE_min2counts[["lux_only_proteins_removep"]]
matched_ALL_DGE_min2counts_VTSUJII <-matched_transcripts_ALL_DGE_min2counts[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_ALL_DGE_min2counts_ALL_transcripts_with_predictedORFs <-matched_transcripts_ALL_DGE_min2counts[["ALL_transcripts_with_predictedORFs"]]

