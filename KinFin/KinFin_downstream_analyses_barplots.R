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
  theme(legend.text = element_text(size = 16))  +           # Legend text
  theme(legend.title = element_text(size = 20))  


data_all_transcripts_barplot +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))



### BCN and Gut Network ###


group_BCN <- c(rep("BCN", 5), rep("Gut" , 5))
origin_of_transcripts_BCN <- rep(c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" ) , 2)
expressed_transcripts_BCN <- c(130, 100, 35, 28, 53, 69, 22, 9, 1, 4) 
data_BCN <- data.frame(group_BCN,origin_of_transcripts_BCN,expressed_transcripts_BCN)


# stacked barplot
data_BCN_barplot<- ggplot(data_BCN, aes(fill=factor(origin_of_transcripts_BCN, levels = c("All Arthropods", "Ostracoda", "Luminini","Luxorina","V.tsujii" )), y=expressed_transcripts_BCN, x=group_BCN)) + 
  geom_bar(position="fill", stat="identity", color="black", size=0.80, width = .60) + 
  scale_fill_manual(values=c('#acacac','#E69F00', '#0072b2','#56B4E9','#CC79A7')) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Number of Expressed Transcripts", x = "Module") +
  labs(fill = "Origin of Genes") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 20))   +
  theme(legend.text = element_text(size = 16))  +           # Legend text
  theme(legend.title = element_text(size = 20))  


data_BCN_barplot  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))



### correct for multiple testing ###


data <- data.frame(
  Group = c("Bio_Upper_Lip", "Compound_eye", "Gut", "DGE"),
  Luminini = c(22, 6, 56, 351),
  ORFs = c(267, 376, 1529, 11669)
)

# perform Fisher's exact tests or chi-square tests for pairwise comparisons
# Bio Upper Lip vs DGE
if (min(table(data[data$Group == "Bio_Upper_Lip", "Luminini"]), table(data[data$Group == "DGE", "Luminini"])) < 5 ||
    min(table(data[data$Group == "Bio_Upper_Lip", "ORFs"]), table(data[data$Group == "DGE", "ORFs"])) < 5) {
  test1 <- fisher.test(matrix(c(data[data$Group == "Bio_Upper_Lip", "Luminini"],
                                data[data$Group == "Bio_Upper_Lip", "ORFs"],
                                data[data$Group == "DGE", "Luminini"],
                                data[data$Group == "DGE", "ORFs"]), nrow = 2))
} else {
  test1 <- chisq.test(matrix(c(data[data$Group == "Bio_Upper_Lip", "Luminini"],
                               data[data$Group == "Bio_Upper_Lip", "ORFs"],
                               data[data$Group == "DGE", "Luminini"],
                               data[data$Group == "DGE", "ORFs"]), nrow = 2))
}

# Compound Eye vs DGE
if (min(table(data[data$Group == "Compound_eye", "Luminini"]), table(data[data$Group == "DGE", "Luminini"])) < 5 ||
    min(table(data[data$Group == "Compound_eye", "ORFs"]), table(data[data$Group == "DGE", "ORFs"])) < 5) {
  test2 <- fisher.test(matrix(c(data[data$Group == "Compound_eye", "Luminini"],
                                data[data$Group == "Compound_eye", "ORFs"],
                                data[data$Group == "DGE", "Luminini"],
                                data[data$Group == "DGE", "ORFs"]), nrow = 2))
} else {
  test2 <- chisq.test(matrix(c(data[data$Group == "Compound_eye", "Luminini"],
                               data[data$Group == "Compound_eye", "ORFs"],
                               data[data$Group == "DGE", "Luminini"],
                               data[data$Group == "DGE", "ORFs"]), nrow = 2))
}

# Gut vs DGE
if (min(table(data[data$Group == "Gut", "Luminini"]), table(data[data$Group == "DGE", "Luminini"])) < 5 ||
    min(table(data[data$Group == "Gut", "ORFs"]), table(data[data$Group == "DGE", "ORFs"])) < 5) {
  test3 <- fisher.test(matrix(c(data[data$Group == "Gut", "Luminini"],
                                data[data$Group == "Gut", "ORFs"],
                                data[data$Group == "DGE", "Luminini"],
                                data[data$Group == "DGE", "ORFs"]), nrow = 2))
} else {
  test3 <- chisq.test(matrix(c(data[data$Group == "Gut", "Luminini"],
                               data[data$Group == "Gut", "ORFs"],
                               data[data$Group == "DGE", "Luminini"],
                               data[data$Group == "DGE", "ORFs"]), nrow = 2))
}

# Correct for multiple testing using Bonferroni correction
p_values <- c(test1$p.value, test2$p.value, test3$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# create a data frame with original and adjusted p-values
results <- data.frame(
  Comparison = c("Bio Upper Lip vs DGE", "Compound Eye vs DGE", "Gut vs DGE"),
  p_value = c(test1$p.value, test2$p.value, test3$p.value),
  adjusted_p_value = adjusted_p_values
)

# view the results
results


