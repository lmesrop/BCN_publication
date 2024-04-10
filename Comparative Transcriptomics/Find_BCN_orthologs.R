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

#import the Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp txt file from the Orthofinder output. 
vargulatsujii_v_skogsbergia_ALL_orthologs <- read.csv("Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp.csv", header = TRUE)

##### subset Vargula_tsujii_cdhit_95.fasta.transdecoder__v__Skogsbergia_sp txt to find all one-to-one orthologs #### 

#take a peak 
head(vargulatsujii_v_skogsbergia_ALL_orthologs)

#first, remove any rows that have more than one ortholog in the V.tsujii column 
vargulatsujii_v_skogsbergia_ALL_orthologs_v1 <- vargulatsujii_v_skogsbergia_ALL_orthologs %>% filter(!grepl(',', Vargula_tsujii_cdhit_95.fasta.transdecoder))

#second, remove any rows that have more than one ortholog in the Skogsbergia sp. column
vargulatsujii_v_skogsbergia_ALL_orthologs_v2 <- vargulatsujii_v_skogsbergia_ALL_orthologs_v1 %>% filter(!grepl(',', Skogsbergia_sp))

#
