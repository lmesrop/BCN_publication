#Description: Barplots for downstream analyses of KinFin output for Vargula tsujii. 
#Author: Lisa Yeter Mesrop 
#Goal: Scripts for barplot figures. 

#load libraries 
library(tidyverse) 
library(matrixStats)
library(dplyr)
library(data.table)
library(ggraph) 
library(graphlayouts)
library(igraph)
library(readxl)
library(plyr)


##### barplots for Figure 4 ######

### Vargula tsujiii ###

group_DE <- c( rep("Light Organ" , 5) , rep("Compound Eye" , 5), rep("Gut" , 5))
expressed_transcripts_DE <- c( 89, 103, 22, 27, 26, 227,67, 6, 20, 56, 983, 249, 56, 27, 214)
origin_of_transcripts_DE <- rep(c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" ) , 3)
data_DE <- data.frame(group_DE,origin_of_transcripts_DE,expressed_transcripts_DE)
data_DE$group_DE <- factor(data_DE$group_DE, levels = c("Light Organ", "Compound Eye", "Gut"))


# stacked barplot
data_DE_barplot<- ggplot(data_DE, aes(fill=factor(origin_of_transcripts_DE, levels = c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" )), y=expressed_transcripts_DE, x= group_DE)) + 
  geom_bar(position="fill", stat="identity", color="black", size=.80) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values=c( '#acacac', '#E69F00', '#0072b2', '#56B4E9', '#CC79A7')) + #'#009E73' 
  labs(y = "Number of Significantly Upregulated Transcripts", x = "Tissue Type") +
  labs(fill = "Origin of Genes") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))   +
  theme(legend.text = element_text(size = 16))  +           # Legend text
  theme(legend.title = element_text(size = 20))  

data_DE_barplot +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))


### All expressed transcripts in DGE Dataset ###

group_all <- c(rep("All Expressed Transcripts" , 5))
origin_of_transcripts_all <- rep(c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" ) , 1)
expressed_transcripts_all <- c(6455, 1691, 351, 358, 2814)
data_all <- data.frame(group_all,origin_of_transcripts_all,expressed_transcripts_all)

# stacked barplot
data_all_transcripts_barplot <- ggplot(data_all, aes(fill=factor(origin_of_transcripts_all, levels = c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" )), y=expressed_transcripts_all, x=group_all)) + 
  geom_bar(position="fill", stat="identity", color="black", size=0.80, width = .2) +
  scale_fill_manual(values=c('#acacac','#E69F00', '#0072b2','#56B4E9','#CC79A7')) +
  scale_y_continuous(labels = scales::percent_format()) +labs(y = "Number of Expressed Transcripts", x = "All Expressed Transcripts") +
  labs(fill = "Origin of Genes")  +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))   +
  theme(legend.text = element_text(size = 16))  +         
  theme(legend.title = element_text(size = 20))  


data_all_transcripts_barplot +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))




