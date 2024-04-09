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

## remove .p1 or .p2 
Vargula_tsujii_cdhit_95_proteins_removep <- Vargula_tsujii_cdhit_95_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
ostra_only_proteins_removep <- ostra_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
arthro_only_proteins_removep <-  arthro_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lum_only_proteins_removep <- lum_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
lux_only_proteins_removep <- lux_only_proteins %>% separate(V1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))

## import dataframes for significantly upregulated genes of the luminous upper lip, compound eye and gut

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

##extract just the names and rename for clarity 

#luminous upper lip 
sigfig_up_bio_upper_lip <-subset(sigOE.annotated_upregulated_logfold_upper_lip_vseye_100122, select = "gene")

#compound eye 
sigfig_up_compound_eye <- subset(sigOE.annotated_downregulated_logfold_upper_lip_vseye_100122, select = "gene") 

#gut 
sigfig_up_gut <- subset(sigOE.annotated_upregulated_logfold_gut_vs_eye_ordered_100122, select = "gene") 

#BCN
BCN_module <- subset(X050222_red_module_gene_annot_power_8_prefilter_5_3, select = "transcript_id") 

#gut module 
gut_module <- subset(vtsujii_purple_gut_module_gene_annot_power_8_prefilter_5_3, select = "transcript_id") 

## create a dataframe for the total number of transcripts that have predicted ORFs (from the OrthoFinder and KinFin analyses).Take the txt files from the KinFin output - arthro_proteins.txt and Vargula_tsujii_cdhit_95_proteins.txt - and combine them into a single dataframe ALL_arthro_and_vtsujii_specific_proteins and remove nested transcripts that have been counted twice and removed .p* 
ALL_arthro_and_vtsujii_specific_proteins_rm <- ALL_arthro_and_vtsujii_specific_proteins %>% distinct()
ALL_arthro_and_vtsujii_specific_proteins_rm_removep <- ALL_arthro_and_vtsujii_specific_proteins_rm %>% separate(...1, c( "transcript_id", "col2"), ".p") %>% select(-c(col2))
#rename for clarity 
ALL_transcripts_with_predictedORFs<- ALL_arthro_and_vtsujii_specific_proteins_rm_removep

#all transcripts with predicted ORFs
ALL_transcripts_with_predictedORFs

#all expressed DGE transcripts (prefilter min 2 counts in more than 5 biological replicates or more). Imported from Juypter Notebook. 
ALL_DGE_min2counts <- subset(dds_count_table_organ_level_prefiltered_strict_assay_ids, select = "transript_ids")

#### function for finding taxonomic restricted genes for each dataset of interest (e.g. significantly upregulated genes and co-expressed genes)

## first, combine the dataframes from the KinFin output for each taxonomic level (w/o nested proteins)

list1 <- list(arthro_only_proteins_removep, ostra_only_proteins_removep, lum_only_proteins_removep, lux_only_proteins_removep, Vargula_tsujii_cdhit_95_proteins_removep, ALL_transcripts_with_predictedORFs)

names_list_1 <- c("arthro_only_proteins_removep", "ostra_only_proteins_removep", "lum_only_proteins_removep", "lux_only_proteins_removep", "Vargula_tsujii_cdhit_95_proteins_removep", "ALL_transcripts_with_predictedORFs")
  
## rename the lists to the corresponding taxonomic level 
names(list1) <- names_list_1

#function to input the dataset of interest and determine how many are restricted to each taxonomic level of interest and print out the number of transcripts for each taxonmic level 

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

## sigfig_up_bio_upper_lip dataset 

## call the function 
matched_transcripts_sigfig_up_bio_upper_lip <- match_transcripts_to_taxonomic_origin(sigfig_up_bio_upper_lip, list1)

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

## sigfig_up_compound_eye dataset 

## call the function 
matched_transcripts_sigfig_up_compound_eye <- match_transcripts_to_taxonomic_origin(sigfig_up_compound_eye, list1)

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

## sigfig_up_gut dataset 

## call the function 
matched_transcripts_sigfig_up_gut <- match_transcripts_to_taxonomic_origin(sigfig_up_gut, list1)

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

# BCN_module dataset 

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

# gut_module dataset 

## call the function 
matched_transcripts_gut_module <- match_transcripts_to_taxonomic_origin(gut_module, list1)

## print the output 

for (df_name in names(matched_transcripts_gut_module)) {
  df <- matched_transcripts_gut_module[[df_name]]
  cat("Dataset: gut_module \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_gut_module_ARTHRO <-matched_transcripts_gut_module[["arthro_only_proteins_removep"]]
matched_gut_module_OSTRA <-matched_transcripts_gut_module[["ostra_only_proteins_removep"]]
matched_gut_module_LUM <- matched_transcripts_gut_module[["lum_only_proteins_removep"]]
matched_gut_module_LUX <-matched_transcripts_gut_module[["lux_only_proteins_removep"]]
matched_gut_module_VTSUJII <-matched_transcripts_gut_module[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_gut_module_ALL_transcripts_with_predictedORFs <-matched_transcripts_gut_module[["ALL_transcripts_with_predictedORFs"]]

# ALL_DGE_min2counts

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

#### find the number of uniquely expressed transcripts at each taxonomic level #### 

## dataframes for the significantly upregulated genes in each tissue
sigfig_up_bio_upper_lip
sigfig_up_compound_eye
sigfig_up_gut

#generate a venn diagram to visualize shared significantly upregulated genes across tissue types
unique_venn_list <- list(
  Bio_Upper_Lip =  sigfig_up_bio_upper_lip$gene , 
  Gut = sigfig_up_gut$gene,
  Compound_Eye = sigfig_up_compound_eye$gene
)

ggvenn(
  unique_venn_list, 
  fill_color = c("#3DB3C2", "#71C23D", "#8E3DC2"),
  stroke_size = .7, set_name_size = 4.5
)

##find the genes uniquely upregulated in the bio upper lip and not in the gut
lightorgan_unique_de  <- sigfig_up_bio_upper_lip$gene[!sigfig_up_bio_upper_lip$gene %in% sigfig_up_gut$gene]
lightorgan_unique_de_df <- as.data.frame(lightorgan_unique_de)
nrow(lightorgan_unique_de_df) 

#unique_de_lightorgan_annot <- Vtsujii_Trinotate_lym_subset %>% filter(transcript_id %in% lightorgan_unique_de_df$lightorgan_unique_de)
#write.csv(unique_de_lightorgan_annot, file = "unique_de_lightorgan_annot.csv" )

##find the genes uniquely upregulated in the gut by subtracting from the bio upper lip 
gut_unique_de_v1 <- sigfig_up_gut$gene[!sigfig_up_gut$gene %in% sigfig_up_bio_upper_lip$gene]
gut_unique_de_v1_df <- as.data.frame(gut_unique_de_v1)
nrow(gut_unique_de_v1_df) 

##now subtract the compound eye from the gut_unique_de_v1_df
gut_unique_de <- gut_unique_de_v1_df$gut_unique_de_v1[!gut_unique_de_v1_df$gut_unique_de_v1 %in% sigfig_up_compound_eye$gene]
gut_unique_de_df <- as.data.frame(gut_unique_de)
nrow(gut_unique_de_df) 

#unique_de_GUT_annot <- Vtsujii_Trinotate_lym_subset %>% filter(transcript_id %in% gut_unique_de_df$gut_unique_de)
#write.csv(unique_de_GUT_annot, file = "unique_de_GUT_annot.csv")

##find the genes uniquely upregulated in the compound eye by subtracting from the gut 
eye_unique_de <- sigfig_up_compound_eye$gene[!sigfig_up_compound_eye$gene %in% sigfig_up_gut$gene]
eye_unique_de_df <- as.data.frame(eye_unique_de)
nrow(eye_unique_de_df) 

#unique_de_EYE_annot <- Vtsujii_Trinotate_lym_subset %>% filter(transcript_id %in% eye_unique_de_df$eye_unique_de)
#write.csv(unique_de_EYE_annot, file = "unique_de_EYE_annot.csv")

## gut_unique_de_df

## call the function 
matched_transcripts_gut_unique_de_df <- match_transcripts_to_taxonomic_origin(gut_unique_de_df, list1)

## print the output 

for (df_name in names(matched_transcripts_gut_unique_de_df)) {
  df <- matched_transcripts_gut_unique_de_df[[df_name]]
  cat("Dataset: gut_unique_de_df \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_ALL_gut_unique_de_df_ARTHRO <-matched_transcripts_gut_unique_de_df[["arthro_only_proteins_removep"]]
matched_ALL_gut_unique_de_df_OSTRA <-matched_transcripts_gut_unique_de_df[["ostra_only_proteins_removep"]]
matched_ALL_gut_unique_de_df_LUM <- matched_transcripts_gut_unique_de_df[["lum_only_proteins_removep"]]
matched_ALL_gut_unique_de_df_LUX <-matched_transcripts_gut_unique_de_df[["lux_only_proteins_removep"]]
matched_ALL_gut_unique_de_df_VTSUJII <-matched_transcripts_gut_unique_de_df[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_ALL_gut_unique_de_df_ALL_transcripts_with_predictedORFs <-matched_transcripts_gut_unique_de_df[["ALL_transcripts_with_predictedORFs"]]

## lightorgan_unique_de_df

## call the function 
matched_transcripts_lightorgan_unique_de_df <- match_transcripts_to_taxonomic_origin(lightorgan_unique_de_df, list1)

## print the output 

for (df_name in names(matched_transcripts_lightorgan_unique_de_df)) {
  df <- matched_transcripts_lightorgan_unique_de_df[[df_name]]
  cat("Dataset: lightorgan_unique_de_df \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_ALL_lightorgan_unique_de_df_ARTHRO <-matched_transcripts_lightorgan_unique_de_df[["arthro_only_proteins_removep"]]
matched_ALL_lightorgan_unique_de_df_OSTRA <-matched_transcripts_lightorgan_unique_de_df[["ostra_only_proteins_removep"]]
matched_ALL_lightorgan_unique_de_df_LUM <- matched_transcripts_lightorgan_unique_de_df[["lum_only_proteins_removep"]]
matched_ALL_lightorgan_unique_de_df_LUX <-matched_transcripts_lightorgan_unique_de_df[["lux_only_proteins_removep"]]
matched_ALL_lightorgan_unique_de_df_VTSUJII <-matched_transcripts_lightorgan_unique_de_df[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_ALL_lightorgan_unique_de_df_ALL_transcripts_with_predictedORFs <-matched_transcripts_lightorgan_unique_de_df[["ALL_transcripts_with_predictedORFs"]]

## eye_unique_de_df

## call the function 
matched_transcripts_eye_unique_de_df <- match_transcripts_to_taxonomic_origin(eye_unique_de_df, list1)

## print the output 

for (df_name in names(matched_transcripts_eye_unique_de_df)) {
  df <- matched_transcripts_eye_unique_de_df[[df_name]]
  cat("Dataset: eye_unique_de_df \n")
  cat("Taxonomic level:", df_name, "\n")
  cat("Number of transcripts:", nrow(df), "\n\n")
}

## save the transcript ids for each taxonomic level 
matched_ALL_eye_unique_de_df_ARTHRO <-matched_transcripts_eye_unique_de_df[["arthro_only_proteins_removep"]]
matched_ALL_eye_unique_de_df_df_OSTRA <-matched_transcripts_eye_unique_de_df[["ostra_only_proteins_removep"]]
matched_ALL_eye_unique_de_df_df_LUM <- matched_transcripts_eye_unique_de_df[["lum_only_proteins_removep"]]
matched_ALL_eye_unique_de_df_LUX <-matched_transcripts_eye_unique_de_df[["lux_only_proteins_removep"]]
matched_ALL_eye_unique_de_df_VTSUJII <-matched_transcripts_eye_unique_de_df[["Vargula_tsujii_cdhit_95_proteins_removep"]]
matched_ALL_eye_unique_de_df_ALL_transcripts_with_predictedORFs <-matched_transcripts_eye_unique_de_df[["ALL_transcripts_with_predictedORFs"]]


