#Description:
#Author:
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


#import the orthologs between V.tsujii and Skogsbergia sp. from the Orthofinder output
vargulatsujii_v_skogsbergia_ALL_orthologs <- read.csv("Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp.csv", header = TRUE)

#take a quick peek 
head(vargulatsujii_v_skogsbergia_ALL_orthologs)

#find all orthogroups that have one-to-one orthologs across the two columns (i.e. first column V.tsujii orthologs and second column Skogsbergia sp.)  
#remove any rows that have multiple orthologs in each column 

#first start with V.tsujii 
vargulatsujii_v_skogsbergia_ALL_orthologs_v1 <- vargulatsujii_v_skogsbergia_ALL_orthologs %>% filter(!grepl(',', Vargula_tsujii_cdhit_95.fasta.transdecoder))

#do the same for Skogsbergia sp. 
vargulatsujii_v_skogsbergia_ALL_orthologs_v2 <- vargulatsujii_v_skogsbergia_ALL_orthologs_v1 %>% filter(!grepl(',', Skogsbergia_sp))

#take a quick peek at the df which now has expression counts of one-to-one orthologs across V.tsujii and Skogsbergia sp. transcriptomes 
head(All_one_to_one_orthlogs_vtsujii_skogs)

#fix the df All_one_to_one_orthlogs_vtsujii_skogs

#remove extra column
All_one_to_one_orthlogs_vtsujii_skogs$X <- NULL
#remove .p1 and .p2 at the end of each transcript ID. TIP: To avoid this issue of multiple open reading frames for one transcript 



