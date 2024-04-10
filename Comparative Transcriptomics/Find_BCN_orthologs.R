#Description: Find BCN one-to-one orthologs and all orthologs across V.tsujii and Skogsbergia sp. transcriptomes. 
#Author: Lisa Yeter Mesrop 
#Goal: Use Orthofinder output to find all BCN one-to-one orthologs and all orthologs. 

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
