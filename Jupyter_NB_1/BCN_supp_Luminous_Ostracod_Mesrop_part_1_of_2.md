<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#-Load-libraries" data-toc-modified-id="-Load-libraries-1"><span class="toc-item-num">1&nbsp;&nbsp;</span> Load libraries</a></span></li><li><span><a href="#-Load-data-" data-toc-modified-id="-Load-data--2"><span class="toc-item-num">2&nbsp;&nbsp;</span> Load data </a></span></li><li><span><a href="#-QC-for-downstream-analyses-" data-toc-modified-id="-QC-for-downstream-analyses--3"><span class="toc-item-num">3&nbsp;&nbsp;</span> QC for downstream analyses </a></span><ul class="toc-item"><li><span><a href="#-Barplot-of-counts-for-all-samples-" data-toc-modified-id="-Barplot-of-counts-for-all-samples--3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span> Barplot of counts for all samples </a></span></li><li><span><a href="#-Remove-samples-with-low-counts-for-WGCNA-" data-toc-modified-id="-Remove-samples-with-low-counts-for-WGCNA--3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span> Remove samples with low counts for WGCNA </a></span></li><li><span><a href="#-Perform-variance-stabilizing-transformation-for-WGCNA-" data-toc-modified-id="-Perform-variance-stabilizing-transformation-for-WGCNA--3.3"><span class="toc-item-num">3.3&nbsp;&nbsp;</span> Perform variance-stabilizing transformation for WGCNA </a></span></li><li><span><a href="#-Visualize-transformed-matrix-with-hierarchical-clustering-and-PCA-" data-toc-modified-id="-Visualize-transformed-matrix-with-hierarchical-clustering-and-PCA--3.4"><span class="toc-item-num">3.4&nbsp;&nbsp;</span> Visualize transformed matrix with hierarchical clustering and PCA </a></span></li></ul></li><li><span><a href="#-WGCNA-" data-toc-modified-id="-WGCNA--4"><span class="toc-item-num">4&nbsp;&nbsp;</span> WGCNA </a></span><ul class="toc-item"><li><span><a href="#-Load-data-" data-toc-modified-id="-Load-data--4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span> Load data </a></span></li><li><span><a href="#-Pick-a-soft-threshold-power-" data-toc-modified-id="-Pick-a-soft-threshold-power--4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span> Pick a soft-threshold power </a></span></li><li><span><a href="#-Run-co-expression-similarity-and-adjacency-and-TOM-" data-toc-modified-id="-Run-co-expression-similarity-and-adjacency-and-TOM--4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span> Run co-expression similarity and adjacency and TOM </a></span></li><li><span><a href="#-Cluster-genes-by-TOM-and-merge-modules-with-very-similar-expression-profiles-" data-toc-modified-id="-Cluster-genes-by-TOM-and-merge-modules-with-very-similar-expression-profiles--4.4"><span class="toc-item-num">4.4&nbsp;&nbsp;</span> Cluster genes by TOM and merge modules with very similar expression profiles </a></span></li><li><span><a href="#-Identify-the-module-that-contains--Vargula-tsujii--c-luciferase
-" data-toc-modified-id="-Identify-the-module-that-contains--Vargula-tsujii--c-luciferase
--4.5"><span class="toc-item-num">4.5&nbsp;&nbsp;</span> Identify the module that contains <em> Vargula tsujii </em> c-luciferase
 </a></span></li><li><span><a href="#-Import-the-annotation-for--Vargula-tsujii--transcriptome-
-" data-toc-modified-id="-Import-the-annotation-for--Vargula-tsujii--transcriptome-
--4.6"><span class="toc-item-num">4.6&nbsp;&nbsp;</span> Import the annotation for <em> Vargula tsujii </em> transcriptome 
 </a></span></li><li><span><a href="#-Identify-other-modules-of-interest-and-incorporate-annotation-information
-" data-toc-modified-id="-Identify-other-modules-of-interest-and-incorporate-annotation-information
--4.7"><span class="toc-item-num">4.7&nbsp;&nbsp;</span> Identify other modules of interest and incorporate annotation information
 </a></span></li><li><span><a href="#-Module-to-trait-heatmap-" data-toc-modified-id="-Module-to-trait-heatmap--4.8"><span class="toc-item-num">4.8&nbsp;&nbsp;</span> Module to trait heatmap </a></span></li></ul></li><li><span><a href="#-Network-visualization-" data-toc-modified-id="-Network-visualization--5"><span class="toc-item-num">5&nbsp;&nbsp;</span> Network visualization </a></span><ul class="toc-item"><li><span><a href="#-Export-BCN-for-Cytoscape--" data-toc-modified-id="-Export-BCN-for-Cytoscape---5.1"><span class="toc-item-num">5.1&nbsp;&nbsp;</span> Export BCN for Cytoscape  </a></span></li></ul></li><li><span><a href="#-GO-enrichment-analyses-for-BCN-" data-toc-modified-id="-GO-enrichment-analyses-for-BCN--6"><span class="toc-item-num">6&nbsp;&nbsp;</span> GO enrichment analyses for BCN </a></span><ul class="toc-item"><li><span><a href="#-Import-GO-annotations-" data-toc-modified-id="-Import-GO-annotations--6.1"><span class="toc-item-num">6.1&nbsp;&nbsp;</span> Import GO annotations </a></span></li><li><span><a href="#-Use-topGO-to-identify-enriched-biological-processes-in-the-BCN--" data-toc-modified-id="-Use-topGO-to-identify-enriched-biological-processes-in-the-BCN---6.2"><span class="toc-item-num">6.2&nbsp;&nbsp;</span> Use topGO to identify enriched biological processes in the BCN  </a></span></li><li><span><a href="#-Use-topGO-to-identify-enriched-molecular-functions-in-BCN--" data-toc-modified-id="-Use-topGO-to-identify-enriched-molecular-functions-in-BCN---6.3"><span class="toc-item-num">6.3&nbsp;&nbsp;</span> Use topGO to identify enriched molecular functions in BCN  </a></span></li></ul></li><li><span><a href="#-BCN-network-connectivity-" data-toc-modified-id="-BCN-network-connectivity--7"><span class="toc-item-num">7&nbsp;&nbsp;</span> BCN network connectivity </a></span><ul class="toc-item"><li><span><a href="#Determine-connectivity-between-intramodular-connectivity-and-module-membership---" data-toc-modified-id="Determine-connectivity-between-intramodular-connectivity-and-module-membership----7.1"><span class="toc-item-num">7.1&nbsp;&nbsp;</span>Determine connectivity between intramodular connectivity and module membership   </a></span></li><li><span><a href="#Identify-genes-in-the-BCN-with-highest-intramodular-connectivity-" data-toc-modified-id="Identify-genes-in-the-BCN-with-highest-intramodular-connectivity--7.2"><span class="toc-item-num">7.2&nbsp;&nbsp;</span>Identify genes in the BCN with highest intramodular connectivity </a></span></li><li><span><a href="#Identify-genes-in-the-bioluminescent-upper-lip-with-highest-intramodular-connectivity-" data-toc-modified-id="Identify-genes-in-the-bioluminescent-upper-lip-with-highest-intramodular-connectivity--7.3"><span class="toc-item-num">7.3&nbsp;&nbsp;</span>Identify genes in the bioluminescent upper lip with highest intramodular connectivity </a></span></li></ul></li><li><span><a href="#-QC-for-Differential-Gene-Expression--" data-toc-modified-id="-QC-for-Differential-Gene-Expression---8"><span class="toc-item-num">8&nbsp;&nbsp;</span> QC for Differential Gene Expression  </a></span></li><li><span><a href="#-Differential-gene-expression-" data-toc-modified-id="-Differential-gene-expression--9"><span class="toc-item-num">9&nbsp;&nbsp;</span> Differential gene expression </a></span><ul class="toc-item"><li><span><a href="#-DGE---Gut-vs-Compound-Eye--" data-toc-modified-id="-DGE---Gut-vs-Compound-Eye---9.1"><span class="toc-item-num">9.1&nbsp;&nbsp;</span> DGE - Gut vs Compound Eye  </a></span></li><li><span><a href="#-DGE---Bioluminescent-Upper-Lip-vs-Compound-Eye--" data-toc-modified-id="-DGE---Bioluminescent-Upper-Lip-vs-Compound-Eye---9.2"><span class="toc-item-num">9.2&nbsp;&nbsp;</span> DGE - Bioluminescent Upper Lip vs Compound Eye  </a></span></li><li><span><a href="#-DGE---Bioluminescent-Upper-Lip-vs-Gut-" data-toc-modified-id="-DGE---Bioluminescent-Upper-Lip-vs-Gut--9.3"><span class="toc-item-num">9.3&nbsp;&nbsp;</span> DGE - Bioluminescent Upper Lip vs Gut </a></span></li><li><span><a href="#-Determine-tissue-specific-expression--" data-toc-modified-id="-Determine-tissue-specific-expression---9.4"><span class="toc-item-num">9.4&nbsp;&nbsp;</span> Determine tissue-specific expression  </a></span></li></ul></li><li><span><a href="#-DGE---GO-enrichment-analyses-for-bioluminescent-upper-lip--" data-toc-modified-id="-DGE---GO-enrichment-analyses-for-bioluminescent-upper-lip---10"><span class="toc-item-num">10&nbsp;&nbsp;</span> DGE - GO enrichment analyses for bioluminescent upper lip  </a></span><ul class="toc-item"><li><span><a href="#-Use-topGO-to-identify-enriched-biological-processes-in-the-bioluminescent-upper-lip--" data-toc-modified-id="-Use-topGO-to-identify-enriched-biological-processes-in-the-bioluminescent-upper-lip---10.1"><span class="toc-item-num">10.1&nbsp;&nbsp;</span> Use topGO to identify enriched biological processes in the bioluminescent upper lip  </a></span></li><li><span><a href="#-Use-topGO-to-identify-enriched-molecular-functions-in-the-bioluminescent-upper-lip-" data-toc-modified-id="-Use-topGO-to-identify-enriched-molecular-functions-in-the-bioluminescent-upper-lip--10.2"><span class="toc-item-num">10.2&nbsp;&nbsp;</span> Use topGO to identify enriched molecular functions in the bioluminescent upper lip </a></span></li></ul></li><li><span><a href="#-References--" data-toc-modified-id="-References---11"><span class="toc-item-num">11&nbsp;&nbsp;</span> References  </a></span></li></ul></div>

**<span style="font-size:30px;"><span style="color:#3A9BDC;"> Supp Material for Luminous Ostracod <em> Vargula tsujii </em> </span>**


The R Juypter Notebook serves as a comprehensive repository encompasing most of the scripts and figures relevant to the bioluminescent ostracod *V.tsujii* analyses outlined in the publication. This notebook includes the following analyses: QC steps, WGCNA, Network Visualization, Network Connectivity, Differential Gene Expression and GO enrichments. 

Author: Lisa Yeter Mesrop 

 

<h1 style="color: #3A9BDC;"> Load libraries</h1> 


```R
#load libraries 
library(WGCNA)
library(tidyverse) 
library(edgeR)
library(matrixStats)
library(DESeq2)
library(dplyr)
library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(repr)
library(topGO)
library(reshape2)
library(plyr)
library(scales)
library(readxl)
library(repr)
library(RColorBrewer)
library(pheatmap)
library(flashClust)
library(ggvenn)


#always use the following WGCNA functions 
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
allowWGCNAThreads(nThreads = 22)

```

<h1 style="color: #3A9BDC;"> Load data </h1> 

The gene expression matrix consists of 57 samples and 41,486 genes before quality control (QC). 


```R
#read in count matrix
merged_counts <- read.csv("merged_single_network.csv", header = TRUE, row.names=1,
                  stringsAsFactors = FALSE)
```


```R
#create the metadata 
meta <- data.frame(row.names = colnames(merged_counts)) 
```


```R
sample_name = c("upper_lip_1","upper_lip_1","upper_lip_1","upper_lip_1","upper_lip_1","upper_lip_1","upper_lip_1","upper_lip_1",
                                  "upper_lip_1","upper_lip_2", "upper_lip_2", "upper_lip_2", "upper_lip_2", "upper_lip_2", 
                                  "eye", "gut", "eye", "gut", "eye", "gut","eye", "gut","eye", "gut","A", "A", "A", "AIF", "AIF", "AIF",
                                  "AII", "AII", "AII", "AIII", "AIII", "AIII", "AIM", "AIM", "AIM", "AIV", "AIV", "AIV","AV", "AV", "AV", "C","C","C","E", "E", "E", "G", "H", "H", "M", "M","M")
```


```R
meta$sample_name <- sample_name 
```


```R
meta$names <- rownames(meta)
```


```R
rownames(meta) <- NULL
```


```R
head(meta)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>sample_name</th><th scope=col>names</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>upper_lip_1</td><td>vtu10.counts.tab</td></tr>
	<tr><th scope=row>2</th><td>upper_lip_1</td><td>vtu1.counts.tab </td></tr>
	<tr><th scope=row>3</th><td>upper_lip_1</td><td>vtu2.counts.tab </td></tr>
	<tr><th scope=row>4</th><td>upper_lip_1</td><td>vtu3.counts.tab </td></tr>
	<tr><th scope=row>5</th><td>upper_lip_1</td><td>vtu4.counts.tab </td></tr>
	<tr><th scope=row>6</th><td>upper_lip_1</td><td>vtu5.counts.tab </td></tr>
</tbody>
</table>



<h1 style="color: #3A9BDC;"> QC for downstream analyses </h1> 


As a preliminary quality control (QC) measure for WGCNA analysis, the overall similarity between samples and transcripts with low counts was assessed, as these counts often introduce noise in co-expression analyses. A filter was applied to the expression matrix, removing transcripts with fewer than 5 counts in more than 3 samples, given that some sample types included a minimum of 3 biological replicates. The QC analyses utilized functions from the DESeq2 package (Love et al., 2014).



```R
#import count table, meta and sample_name objects into a DESeq2 object 
dds_count_table <- DESeqDataSetFromMatrix(countData = merged_counts, colData = meta, design = ~sample_name)
```

<h2 style="color: #3A9BDC;"> Barplot of counts for all samples </h2> 

The counts of filtered reads were examined using a bar plot.


```R
#the number of prefiltered counts for each sample 
colSums(assay(dds_count_table))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>vtu10.counts.tab</dt><dd>881916</dd><dt>vtu1.counts.tab</dt><dd>953139</dd><dt>vtu2.counts.tab</dt><dd>429668</dd><dt>vtu3.counts.tab</dt><dd>882596</dd><dt>vtu4.counts.tab</dt><dd>887784</dd><dt>vtu5.counts.tab</dt><dd>717749</dd><dt>vtu7.counts.tab</dt><dd>1433938</dd><dt>vtu8.counts.tab</dt><dd>488033</dd><dt>vtu9.counts.tab</dt><dd>896901</dd><dt>Vt.1A.counts.tab</dt><dd>1826043</dd><dt>Vt.2A.counts.tab</dt><dd>1388482</dd><dt>Vt.3A.counts.tab</dt><dd>1483095</dd><dt>Vt.4A.counts.tab</dt><dd>1214013</dd><dt>Vt.5A.counts.tab</dt><dd>1636430</dd><dt>Vt.1B.counts.tab</dt><dd>798885</dd><dt>Vt.1C.counts.tab</dt><dd>1365090</dd><dt>Vt.2B.counts.tab</dt><dd>342949</dd><dt>Vt.2C.counts.tab</dt><dd>1927280</dd><dt>Vt.3B.counts.tab</dt><dd>363212</dd><dt>Vt.3C.counts.tab</dt><dd>2470404</dd><dt>Vt.4B.counts.tab</dt><dd>943558</dd><dt>Vt.4C.counts.tab</dt><dd>1623438</dd><dt>Vt.5B.counts.tab</dt><dd>794665</dd><dt>Vt.5C.counts.tab</dt><dd>2014035</dd><dt>A.1.counts.tab</dt><dd>995203</dd><dt>A.2.counts.tab</dt><dd>583498</dd><dt>A.3.counts.tab</dt><dd>683680</dd><dt>AIF.1.counts.tab</dt><dd>1331646</dd><dt>AIF.2.counts.tab</dt><dd>1253052</dd><dt>AIF.3.counts.tab</dt><dd>688206</dd><dt>AII.1.counts.tab</dt><dd>993596</dd><dt>AII.2.counts.tab</dt><dd>808877</dd><dt>AII.3.counts.tab</dt><dd>976551</dd><dt>AIII.1.counts.tab</dt><dd>916076</dd><dt>AIII.2.counts.tab</dt><dd>765084</dd><dt>AIII.3.counts.tab</dt><dd>686605</dd><dt>AIM.1.counts.tab</dt><dd>906177</dd><dt>AIM.2.counts.tab</dt><dd>748870</dd><dt>AIM.3.counts.tab</dt><dd>910231</dd><dt>AIV.1.counts.tab</dt><dd>594556</dd><dt>AIV.2.counts.tab</dt><dd>741113</dd><dt>AIV.3.counts.tab</dt><dd>660913</dd><dt>AV.1.counts.tab</dt><dd>76179</dd><dt>AV.2.counts.tab</dt><dd>58955</dd><dt>AV.4.counts.tab</dt><dd>132446</dd><dt>C.1.counts.tab</dt><dd>969927</dd><dt>C.2.counts.tab</dt><dd>763769</dd><dt>C.3.counts.tab</dt><dd>851898</dd><dt>E.1.counts.tab</dt><dd>1054740</dd><dt>E.2.counts.tab</dt><dd>1084075</dd><dt>E.3.counts.tab</dt><dd>831205</dd><dt>G.1.counts.tab</dt><dd>1210799</dd><dt>H.1.counts.tab</dt><dd>815551</dd><dt>H.2.counts.tab</dt><dd>782567</dd><dt>M1.1.counts.tab</dt><dd>102298</dd><dt>M2.1.counts.tab</dt><dd>118595</dd><dt>M3.1.counts.tab</dt><dd>48696</dd></dl>




```R
#function from the R package repr to visualize figures in Jupyter Notebook(optional)
options(repr.plot.width=10, repr.plot.height=10)
```


```R
#extract sample names for barplot 
names <- dds_count_table$sample_name
```


```R
#visualize prefiltered raw counts in a barplot  

librarySizes <- colSums(assay(dds_count_table))

par(mar=c(10,5,2,2))  
barplot(librarySizes, 
        names=names, 
        las=2, 
        cex.names=.5) +title(main = "Barplot of Count Distributions of Samples", line = -1, outer = TRUE)
       


```






![png](output_18_1.png)


<h2 style="color: #3A9BDC;"> Remove samples with low counts for WGCNA </h2> 

Based on the bar plot count distribution in Section 3.1, the following six samples were removed due to very low counts: M1.1.counts.tab, M2.1.counts.tab, M3.1.counts.tab, AV.1.counts.tab, AV.2.counts.tab, AV.4.counts.tab.



```R
#update the merged counts, meta and sample_name objects with the six samples removed 
merged_counts_subset <- subset(merged_counts, select= -c(M1.1.counts.tab, M2.1.counts.tab, M3.1.counts.tab, AV.1.counts.tab, AV.2.counts.tab, AV.4.counts.tab))
meta_subset <- data.frame(row.names = colnames(merged_counts_subset)) 
sample_name_subset = c("upper_lip","upper_lip","upper_lip","upper_lip","upper_lip","upper_lip","upper_lip","upper_lip",
                                  "upper_lip","upper_lip", "upper_lip", "upper_lip", "upper_lip", "upper_lip", 
                                  "eye", "gut", "eye", "gut", "eye", "gut","eye", "gut","eye", "gut","A", "A", "A", "AIF", "AIF", "AIF",
                                  "AII", "AII", "AII", "AIII", "AIII", "AIII", "AIM", "AIM", "AIM", "AIV", "AIV", "AIV", "C","C","C","E", "E", "E", "G", "H", "H")


meta_subset$sample_name_subset <- sample_name_subset
meta_subset$names_subset <- rownames(meta_subset)
rownames(meta_subset) <- NULL

```


```R
nrow(meta_subset)
```


51


 <h2 style="color: #3A9BDC;"> Perform variance-stabilizing transformation for WGCNA </h2> 
 
Perform variance-stabilizing transformation using DEseq2 package (Love et al., 2014). The authors of WGCNA recommend variance-stabilizing transformation as a normalization method before conducting network analyses (Langfelder and Horvath, 2008).



```R
#import updated count table, meta and sample_name objects into a DESeq2 object 
dds_count_table_subset<- DESeqDataSetFromMatrix(countData = merged_counts_subset, colData = meta_subset, design = ~sample_name_subset)
```


```R
#filter the count table 
dds_merged_table_prefiltered <- dds_count_table_subset[rowSums(counts(dds_count_table_subset) >=5) >=3,]
```


```R
nrow(dds_merged_table_prefiltered)
```


31097



```R
#run DESeq2 function with variance-stabilizing transformation 
dds_subset <- DESeq(dds_merged_table_prefiltered, betaPrior = FALSE, parallel = TRUE)
#perform a variance-stabilizing transformation
vsd_subset <- varianceStabilizingTransformation(dds_subset)


```

    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates: 38 workers
    
    mean-dispersion relationship
    
    final dispersion estimates, fitting model and testing: 38 workers
    
    -- replacing outliers and refitting for 7985 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)
    
    estimating dispersions
    
    fitting model and testing
    


<h2 style="color: #3A9BDC;"> Visualize transformed matrix with hierarchical clustering and PCA </h2> 


```R
#transpose the matrix 
sampleDists_subset <- dist(t(assay(vsd_subset)))

```


```R
#perform heatmap
sampleDistMatrix <- as.matrix(sampleDists_subset)
rownames(sampleDistMatrix) <- paste(colData(dds_merged_table_prefiltered)$sample_name_subset) 
colnames(sampleDistMatrix) <- colData(dds_merged_table_prefiltered)$names_subset
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
          clustering_distance_rows=sampleDists_subset,
         clustering_distance_cols=sampleDists_subset,
         col=colors)

```


![png](output_29_0.png)



```R
#perform PCA 

pcaData <- plotPCA(vsd_subset, intgroup="sample_name_subset", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample_name_subset, shape=sample_name_subset)) +
  geom_point(size=5) + 
labs(color = "Sample Types")+ labs(shape = "Sample Types")+
scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))+
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme(
    panel.grid.major = element_line(colour = "gray97", size = 1),
    panel.grid.minor = element_line(linetype = "dotted"),
    panel.background = element_rect(fill = NA),
    legend.key = element_rect(fill = "gray100"),
    axis.line = element_line(size = 0.5, linetype = "solid"),
   panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + theme(
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16), 
  legend.title = element_text(size = 16), 
  legend.text = element_text(size = 12)    
)


```


![png](output_30_0.png)


<h1 style="color: #3A9BDC;"> WGCNA </h1> 


After completing the prefiltering and normalization steps above, the count matrix now includes 51 samples and 31,097 transcript and is ready for WGCNA. WGCNA (Weighted Gene Co-expression Network Analysis) identifies clusters (modules) of highly correlated genes by constructing a network based on pairwise correlations between gene expression profiles (Langfelder and Horvath, 2008). These modules often correspond to specific biological processes or pathways, indicating that the genes within a module may be part of the same regulatory process. By analyzing the connectivity of genes within each module, WGCNA also helps identify key drivers or hub genes, providing insights into gene regulation and the biological processes they govern. The following scripts are from the WGCNA package (Langfelder and Horvath, 2008). 


<h2 style="color: #3A9BDC;"> Load data </h2> 



```R
#use the prefiltered and normalized count matrix from the DEseq2
vsd_subset_matrix <- assay(vsd_subset)
## many functions expect the matrix to be transposed
datExpr <- t(vsd_subset_matrix) 
## check rows/cols
nrow(datExpr)
ncol(datExpr)
```


51



31097



```R
#ready to start WGCNA Analysis 

#run this to check if there are gene outliers
gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 
```

     Flagging genes and samples with too many missing values...
      ..step 1



TRUE


<h2 style="color: #3A9BDC;"> Pick a soft-threshold power </h2> 



```R
# choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

```

    pickSoftThreshold: will use block size 1438.
     pickSoftThreshold: calculating connectivity for given powers...
       ..working on genes 1 through 1438 of 31097
       ..working on genes 1439 through 2876 of 31097
       ..working on genes 2877 through 4314 of 31097
       ..working on genes 4315 through 5752 of 31097
       ..working on genes 5753 through 7190 of 31097
       ..working on genes 7191 through 8628 of 31097
       ..working on genes 8629 through 10066 of 31097
       ..working on genes 10067 through 11504 of 31097
       ..working on genes 11505 through 12942 of 31097
       ..working on genes 12943 through 14380 of 31097
       ..working on genes 14381 through 15818 of 31097
       ..working on genes 15819 through 17256 of 31097
       ..working on genes 17257 through 18694 of 31097
       ..working on genes 18695 through 20132 of 31097
       ..working on genes 20133 through 21570 of 31097
       ..working on genes 21571 through 23008 of 31097
       ..working on genes 23009 through 24446 of 31097
       ..working on genes 24447 through 25884 of 31097
       ..working on genes 25885 through 27322 of 31097
       ..working on genes 27323 through 28760 of 31097
       ..working on genes 28761 through 30198 of 31097
       ..working on genes 30199 through 31097 of 31097
       Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    1      1    0.157 -0.658          0.818  6860.0  6.55e+03  12600
    2      2    0.787 -1.170          0.915  2420.0  1.96e+03   6880
    3      3    0.894 -1.170          0.885  1110.0  6.97e+02   4360
    4      4    0.929 -1.120          0.912   617.0  2.78e+02   3180
    5      5    0.943 -1.040          0.981   392.0  1.21e+02   2500
    6      6    0.940 -1.010          0.991   275.0  5.58e+01   2150
    7      7    0.943 -0.978          0.994   207.0  2.74e+01   1910
    8      8    0.947 -0.957          0.985   164.0  1.41e+01   1740
    9      9    0.948 -0.943          0.970   135.0  7.49e+00   1610
    10    10    0.951 -0.929          0.960   114.0  4.12e+00   1500
    11    12    0.950 -0.914          0.947    86.6  1.37e+00   1330
    12    14    0.936 -0.916          0.919    68.8  4.98e-01   1200
    13    16    0.922 -0.920          0.900    56.4  1.96e-01   1090
    14    18    0.923 -0.923          0.902    47.3  8.15e-02    999
    15    20    0.914 -0.930          0.891    40.3  3.58e-02    920



```R
# plot the results:
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

```


![png](output_37_0.png)



![png](output_37_1.png)


<h2 style="color: #3A9BDC;"> Run co-expression similarity and adjacency and TOM </h2> 



```R
# co-expression similarity and adjacency using assigned softpower
softPower=8
adjacency = adjacency(datExpr, power = softPower)
```


```R
# calculate the Topological Overlap Matrix (TOM)
# turn adjacency into topological overlap, i.e. translate the adjacency into 
# topological overlap matrix and calculate the corresponding dissimilarity 
TOM = TOMsimilarity(adjacency, TOMType = "signed", verbose = 5);
dissTOM = 1-TOM;
```

    ..connectivity..
    ..matrix multiplication (system BLAS)..
    ..normalization..
    ..done.


<h2 style="color: #3A9BDC;"> Cluster genes by TOM and merge modules with very similar expression profiles </h2> 



```R
# call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

```


```R
# plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


```


![png](output_43_0.png)



```R
# set the min module size
minModuleSize = 30;
# module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
```

     ..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.
     ..done.



```R
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)


```


    dynamicMods
        0     1     2     3     4     5     6     7     8     9    10    11    12 
     6310 10309  7139  2726  1603   944   890   332   289   188   178    82    56 
       13 
       51 



```R
table(dynamicColors)
```


    dynamicColors
          black        blue       brown       green greenyellow        grey 
            332        7139        2726         944          82        6310 
        magenta        pink      purple         red      salmon         tan 
            188         289         178         890          51          56 
      turquoise      yellow 
          10309        1603 



```R
# plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


```


![png](output_47_0.png)



```R
# calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")


```


![png](output_48_0.png)



```R
# call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# merged module colors
mergedColors = merge$colors;
# eigengenes of the new merged modules 
mergedMEs = merge$newMEs;

```

     mergeCloseModules: Merging modules whose distance is less than 0.25
       multiSetMEs: Calculating module MEs.
         Working on set 1 ...
         moduleEigengenes: Calculating 14 module eigengenes in given set.
       multiSetMEs: Calculating module MEs.
         Working on set 1 ...
         moduleEigengenes: Calculating 13 module eigengenes in given set.
       Calculating new MEs...
       multiSetMEs: Calculating module MEs.
         Working on set 1 ...
         moduleEigengenes: Calculating 13 module eigengenes in given set.



```R
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```


![png](output_50_0.png)


<h2 style="color: #3A9BDC;"> Identify the module that contains <em> Vargula tsujii </em> c-luciferase
 </h2> 
 
The functionally demonstrated c-luciferase gene, along with some phylogenetically similar c-luciferase genes, are located in the red module, which we refer to as the Bioluminescent Co-Regulatory Network (BCN). 


```R
#determine which column number in the datExpr object that corresponds to c-luciferase
which(colnames(datExpr) == "NODE_10321_length_1972_cov_1770.41_g3092_i1") 
```


389



```R
#determine the color of the module with c-luciferase
dynamicColors[[389]]
```


'red'



```R
# extract all modules as a table for downstream analyses 
SubGeneNames = colnames(datExpr)
for (color in dynamicColors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

```


```R
# extract the red module of interest from the dynamicColors object 
SubGeneNames = colnames(datExpr)
red=as.data.frame(SubGeneNames[which(dynamicColors=="red")])
names(red)[1] <- "transcript_id"
```


```R
head(red)
```


<table class="dataframe">
<caption>A data.frame: 6 × 1</caption>
<thead>
	<tr><th></th><th scope=col>transcript_id</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_100020_length_140_cov_8.6_g92800_i0   </td></tr>
	<tr><th scope=row>2</th><td>NODE_10035_length_2011_cov_7.68814_g7072_i0</td></tr>
	<tr><th scope=row>3</th><td>NODE_1015_length_4889_cov_22.1146_g714_i0  </td></tr>
	<tr><th scope=row>4</th><td>NODE_10218_length_1985_cov_18.5554_g7198_i0</td></tr>
	<tr><th scope=row>5</th><td>NODE_102839_length_135_cov_19.3_g95619_i0  </td></tr>
	<tr><th scope=row>6</th><td>NODE_10321_length_1972_cov_1770.41_g3092_i1</td></tr>
</tbody>
</table>



<h2 style="color: #3A9BDC;"> Import the annotation for <em> Vargula tsujii </em> transcriptome 
 </h2> 
  


```R
#read in the Trinotate sheet for Vargula tsujii 
Trinotate_lym_subset <- read_excel("Trinotate_lym_subset.xlsx")

```

 <h2 style="color: #3A9BDC;"> Identify other modules of interest and incorporate annotation information
 </h2> 
 



```R
#incorporate annotation for the red module (the BCN)
red_module_gene_annot <- setDT(Trinotate_lym_subset, key = 'transcript_id')[J(red)]

```


```R
# red module (the BCN)
head(red_module_gene_annot)
```


<table class="dataframe">
<caption>A data.table: 6 × 16</caption>
<thead>
	<tr><th scope=col>#gene_id</th><th scope=col>transcript_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>prot_id</th><th scope=col>prot_coords</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>g92800</td><td>NODE_100020_length_140_cov_8.6_g92800_i0   </td><td>.                                                                                                                                                                                                                                                                                                                                                </td><td>.</td><td>.                                             </td><td>.           </td><td>.                                                                                                                                                                                                                                                                                                                                            </td><td>.                                                                                                                                                                                  </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>g7072 </td><td>NODE_10035_length_2011_cov_7.68814_g7072_i0</td><td>.                                                                                                                                                                                                                                                                                                                                                </td><td>.</td><td>.                                             </td><td>.           </td><td>.                                                                                                                                                                                                                                                                                                                                            </td><td>.                                                                                                                                                                                  </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>g714  </td><td>NODE_1015_length_4889_cov_22.1146_g714_i0  </td><td>ROA1_SCHAM^ROA1_SCHAM^Q:4777-4220,H:1-186^67.725%ID^E:7.26e-70^RecName: Full=Heterogeneous nuclear ribonucleoprotein A1, A2/B1 homolog;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Polyneoptera; Orthoptera; Caelifera; Acrididea; Acridomorpha; Acridoidea; Acrididae; Cyrtacanthacridinae; Schistocerca</td><td>.</td><td>NODE_1015_length_4889_cov_22.1146_g714_i0.p1  </td><td>4016-4768[-]</td><td>ROA1_SCHAM^ROA1_SCHAM^Q:1-183,H:4-186^68.817%ID^E:2.58e-86^RecName: Full=Heterogeneous nuclear ribonucleoprotein A1, A2/B1 homolog;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Polyneoptera; Orthoptera; Caelifera; Acrididea; Acridomorpha; Acridoidea; Acrididae; Cyrtacanthacridinae; Schistocerca</td><td>PF00076.22^RRM_1^RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)^14-81^E:7.4e-17`PF00076.22^RRM_1^RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)^105-175^E:1.5e-17</td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>GO:0005634^cellular_component^nucleus`GO:0003729^molecular_function^mRNA binding`GO:0007281^biological_process^germ cell development                                                                                                                                                                                                                                  </td><td>GO:0003676^molecular_function^nucleic acid binding                     </td><td>.</td><td>.</td></tr>
	<tr><td>g7198 </td><td>NODE_10218_length_1985_cov_18.5554_g7198_i0</td><td>WSD_ACIAD^WSD_ACIAD^Q:1368-265,H:62-446^21.76%ID^E:1.15e-12^RecName: Full=O-acyltransferase WSD;^Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter                                                                                                                                                    </td><td>.</td><td>NODE_10218_length_1985_cov_18.5554_g7198_i0.p1</td><td>256-1953[-] </td><td>WSD_ACIAD^WSD_ACIAD^Q:196-563,H:62-446^22.005%ID^E:9.38e-13^RecName: Full=O-acyltransferase WSD;^Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter                                                                                                                                                </td><td>PF06974.13^DUF1298^Protein of unknown function (DUF1298)^418-558^E:2.4e-27                                                                                                         </td><td>.                  </td><td>ExpAA=83.94^PredHel=4^Topology=i36-58o73-95i452-474o504-526i</td><td>ENOG410XS7A^Acyltransferase, ws dgat mgat</td><td>KEGG:aci:ACIAD0832`KO:K00635</td><td>GO:0102966^molecular_function^arachidoyl-CoA:1-dodecanol O-acyltransferase activity`GO:0004144^molecular_function^diacylglycerol O-acyltransferase activity`GO:0047196^molecular_function^long-chain-alcohol O-fatty-acyltransferase activity`GO:0006071^biological_process^glycerol metabolic process`GO:0019432^biological_process^triglyceride biosynthetic process</td><td>GO:0004144^molecular_function^diacylglycerol O-acyltransferase activity</td><td>.</td><td>.</td></tr>
	<tr><td>g95619</td><td>NODE_102839_length_135_cov_19.3_g95619_i0  </td><td>.                                                                                                                                                                                                                                                                                                                                                </td><td>.</td><td>.                                             </td><td>.           </td><td>.                                                                                                                                                                                                                                                                                                                                            </td><td>.                                                                                                                                                                                  </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>g3092 </td><td>NODE_10321_length_1972_cov_1770.41_g3092_i1</td><td>LUCI_VARHI^LUCI_VARHI^Q:106-1749,H:2-555^47.227%ID^E:0^RecName: Full=Luciferin 2-monooxygenase;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Crustacea; Oligostraca; Ostracoda; Myodocopa; Myodocopida; Cypridinoidea; Cypridinidae; Vargula                                                                                                       </td><td>.</td><td>NODE_10321_length_1972_cov_1770.41_g3092_i1.p1</td><td>94-1752[+]  </td><td>LUCI_VARHI^LUCI_VARHI^Q:5-552,H:2-555^47.227%ID^E:0^RecName: Full=Luciferin 2-monooxygenase;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Crustacea; Oligostraca; Ostracoda; Myodocopa; Myodocopida; Cypridinoidea; Cypridinidae; Vargula                                                                                                      </td><td>PF00094.25^VWD^von Willebrand factor type D domain^80-234^E:2.7e-14`PF00094.25^VWD^von Willebrand factor type D domain^319-460^E:3.6e-21                                           </td><td>sigP:1^18^0.906^YES</td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>GO:0016831^molecular_function^carboxy-lyase activity`GO:0047712^molecular_function^Cypridina-luciferin 2-monooxygenase activity`GO:0008218^biological_process^bioluminescence                                                                                                                                                                                         </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#extract the module that is significantly correlated with the gut sample (from Section 4.8) from the dynamicColors object   
purple_gut=as.data.frame(SubGeneNames[which(dynamicColors=="purple")])
names(purple_gut)[1] <- "transcript_id"
purple_gut_gene_annot <- setDT(Trinotate_lym_subset, key = 'transcript_id')[J(purple_gut)]
```


```R
#extract the module that is correlated with forced lum sample (from Section 4.8) from the dynamicColors object   
SubGeneNames = colnames(datExpr)
salmon_forced_lum=as.data.frame(SubGeneNames[which(dynamicColors=="salmon")])
names(salmon_forced_lum)[1] <- "transcript_id"

#match the annotation with each transcript 
salmon_forced_lum_module_gene_annot <- setDT(Trinotate_lym_subset, key = 'transcript_id')[J(salmon_forced_lum)]

```

<h2 style="color: #3A9BDC;"> Module to trait heatmap </h2> 

Identify modules (network) that are significantly associated with samples. The module eigengene provides a representative measure of the gene expression patterns within a module, allowing correlation with these traits to determine the most significant associations (Langfelder and Horvath, 2008). This correlation can then be visualized in a module-to-trait heatmap to determine the most significant associations.The following scripts are from the WGCNA package (Langfelder and Horvath, 2008).  

 


```R
traitData <- read_excel("traits_all_tsujii_updated.xlsx")
```


```R
traitData <- as.data.frame(traitData)
```


```R
names(traitData)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Samples'</li><li>'Samples_names'</li><li>'Bio_UpperLip'</li><li>'Gut'</li><li>'Eye'</li><li>'Instar_Stage_A'</li><li>'Instar_Stage_AI'</li><li>'Instar_Stage_AII'</li><li>'Instar_Stage_AIII'</li><li>'Instar_Stage_AIV'</li><li>'Adult_Morning_Timepoint'</li><li>'Adult_forced_biolum'</li></ol>




```R
head(traitData)
```


<table class="dataframe">
<caption>A data.frame: 6 × 12</caption>
<thead>
	<tr><th></th><th scope=col>Samples</th><th scope=col>Samples_names</th><th scope=col>Bio_UpperLip</th><th scope=col>Gut</th><th scope=col>Eye</th><th scope=col>Instar_Stage_A</th><th scope=col>Instar_Stage_AI</th><th scope=col>Instar_Stage_AII</th><th scope=col>Instar_Stage_AIII</th><th scope=col>Instar_Stage_AIV</th><th scope=col>Adult_Morning_Timepoint</th><th scope=col>Adult_forced_biolum</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>A.1.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>2</th><td>A.2.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>3</th><td>A.3.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>AIF.1.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>5</th><td>AIF.2.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>AIF.3.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
allTraits = traitData
```


```R
head(allTraits)
```


<table class="dataframe">
<caption>A data.frame: 6 × 12</caption>
<thead>
	<tr><th></th><th scope=col>Samples</th><th scope=col>Samples_names</th><th scope=col>Bio_UpperLip</th><th scope=col>Gut</th><th scope=col>Eye</th><th scope=col>Instar_Stage_A</th><th scope=col>Instar_Stage_AI</th><th scope=col>Instar_Stage_AII</th><th scope=col>Instar_Stage_AIII</th><th scope=col>Instar_Stage_AIV</th><th scope=col>Adult_Morning_Timepoint</th><th scope=col>Adult_forced_biolum</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>A.1.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>2</th><td>A.2.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>3</th><td>A.3.counts.tab  </td><td>A  </td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>AIF.1.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>5</th><td>AIF.2.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>AIF.3.counts.tab</td><td>AIF</td><td>0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
names(allTraits)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Samples'</li><li>'Samples_names'</li><li>'Bio_UpperLip'</li><li>'Gut'</li><li>'Eye'</li><li>'Instar_Stage_A'</li><li>'Instar_Stage_AI'</li><li>'Instar_Stage_AII'</li><li>'Instar_Stage_AIII'</li><li>'Instar_Stage_AIV'</li><li>'Adult_Morning_Timepoint'</li><li>'Adult_forced_biolum'</li></ol>




```R
All_samples = rownames(datExpr)
```


```R
traitRows = match(All_samples, allTraits$Samples)
```


```R
traitRows
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>28</li><li>29</li><li>30</li><li>31</li><li>32</li><li>33</li><li>34</li><li>35</li><li>36</li><li>37</li><li>38</li><li>39</li><li>40</li><li>41</li><li>42</li><li>43</li><li>44</li><li>45</li><li>46</li><li>47</li><li>48</li><li>49</li><li>50</li><li>51</li><li>1</li><li>2</li><li>3</li><li>4</li><li>5</li><li>6</li><li>7</li><li>8</li><li>9</li><li>10</li><li>11</li><li>12</li><li>13</li><li>14</li><li>15</li><li>16</li><li>17</li><li>18</li><li>19</li><li>20</li><li>21</li><li>22</li><li>23</li><li>24</li><li>25</li><li>26</li><li>27</li></ol>




```R
datTraits = allTraits[traitRows, -1]
```


```R
rownames(datTraits) = allTraits[traitRows, 1]

```


```R
# double check that row names agree
table(rownames(datTraits)==rownames(datExpr))
```


    
    TRUE 
      51 



```R
datTraits$Samples_names <- NULL
```


```R
head(datTraits)
```


<table class="dataframe">
<caption>A data.frame: 6 × 10</caption>
<thead>
	<tr><th></th><th scope=col>Bio_UpperLip</th><th scope=col>Gut</th><th scope=col>Eye</th><th scope=col>Instar_Stage_A</th><th scope=col>Instar_Stage_AI</th><th scope=col>Instar_Stage_AII</th><th scope=col>Instar_Stage_AIII</th><th scope=col>Instar_Stage_AIV</th><th scope=col>Adult_Morning_Timepoint</th><th scope=col>Adult_forced_biolum</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>vtu10.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>vtu1.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>vtu2.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>vtu3.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>vtu4.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>vtu5.counts.tab</th><td>1</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
# sample network based on squared Euclidean distance
A=adjacency(t(datExpr),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
```


```R
# designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# convert traits to a color representation:
# where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits),"C",sep="")
datColors=data.frame(outlierC=outlierColor,traitColors)
# plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and trait heatmap")
```


![png](output_81_0.png)



```R
# choose a module assignment
moduleColors_Vtsujii=dynamicColors
# define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr,moduleColors_Vtsujii)$eigengenes
MEsFemale = orderMEs(MEs0)
modTraitCor = cor(MEsFemale, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
#color code each association by the correlation value; display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",
signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(8, 8.5, 3, 3))
# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits),
yLabels = names(MEsFemale), ySymbols = names(MEsFemale),
colorLabels =FALSE,colors = blueWhiteRed(50),textMatrix=textMatrix,
setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
main = paste("Module-trait relationships"))
```


![png](output_82_0.png)


<h1 style="color: #3A9BDC;"> Network visualization </h1> 

The Bioluminescent Co-regulatory Network (BCN) has a total of 889 co-expressed transcripts making it challenging to visualize the entire network topology comprehensively. To address this, we used Cytoscape to visualize the network and we subsetted the network to only include genes that are signifcantly upregulated (uniquely) in the bioluminescent upper lip from Section 9.4. 

<h2 style="color: #3A9BDC;"> Export BCN for Cytoscape  </h2> 

Export the network of interest into edge and node list files for Cytoscape (Shannon et al., 2003). 


```R
# select module of interest
module = "red"
# select module probes
genes = colnames(datExpr)
inModule = dynamicColors==module
modProbes = genes[inModule]
# select the corresponding Topological Overlap
modTOM = TOM_visnet_power8[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes) 

# export the network into edge and node list files Cytoscape 
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
altNodeNames = modProbes,
nodeAttr = dynamicColors[inModule]) 
```


```R
#determine how many of the significantly upregulated genes expressed in the bioluminescent upper lip (from Section 9.4) are found in the BCN 
BCN_DE_unique_Bio_UpperLip  <- red  %>%
  filter(transcript_id %in% unique_genes_bio_upper_lip_info_annot$transcript_id)

#add annotation 
BCN_DE_unique_Bio_UpperLip_annot <- left_join(BCN_DE_unique_Bio_UpperLip,Trinotate_lym_subset,by="transcript_id") 

```

 <h1 style="color: #3A9BDC;"> GO enrichment analyses for BCN </h1> 


GO enrichment analyses was performed using topGO package (Alexa Rahnenfuhrer, 2024). The figures representing the GO enrichments presented here were utilized for analysis, but were not included in the main text of the publication. The GO enrichments were visualized with GO-Figure! package for publication (Reijnders and Waterhouse, 2021). GO-Figure! reference https://github.com/lmesrop/BCN_publication/tree/main/Go-Figure!. 


 

 <h2 style="color: #3A9BDC;"> Import GO annotations </h2> 




```R
#import the trinotate go sheet from the Trinotate output 
geneID2GO <- readMappings(file ="Trinotate_go_lym.txt")
geneNames <- names(geneID2GO)

```


```R
#save the transcript ids of all the annotated genes under geneNames object 
geneNames<- as.character(Trinotate_lym_subset$transcript_id)
#save the transcript 
myInterestingGenes= as.character(red[,1])

```


```R
#subset the genesNames by the transcript IDs in my red module 
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
head(geneList)

```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>NODE_100000_length_141_cov_1.05814_g92780_i0</dt><dd>0</dd><dt>NODE_100001_length_141_cov_1.03488_g92781_i0</dt><dd>0</dd><dt>NODE_100002_length_141_cov_0.325581_g92782_i0</dt><dd>0</dd><dt>NODE_100003_length_141_cov_0.0116279_g92783_i0</dt><dd>0</dd><dt>NODE_100004_length_141_cov_0_g92784_i0</dt><dd>0</dd><dt>NODE_100005_length_141_cov_0_g92785_i0</dt><dd>0</dd></dl>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'0'</li><li>'1'</li></ol>
</details>



<h2 style="color: #3A9BDC;"> Use topGO to identify enriched biological processes in the BCN  </h2> 



```R
#run the topGO function. 
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO)

results_go <- runTest(GOdata, algorithm="weight01", statistic="Fisher")


```

    
    Building most specific GOs .....
    
    	( 12645 GO terms found. )
    
    
    Build GO DAG topology ..........
    
    	( 13427 GO terms and 31130 relations. )
    
    
    Annotating nodes ...............
    
    	( 15506 genes annotated to the GO terms. )
    
    
    			 -- Weight01 Algorithm -- 
    
    		 the algorithm is scoring 2181 nontrivial nodes
    		 parameters: 
    			 test statistic: fisher
    
    
    	 Level 16:	1 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 15:	6 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 14:	11 nodes to be scored	(6 eliminated genes)
    
    
    	 Level 13:	30 nodes to be scored	(210 eliminated genes)
    
    
    	 Level 12:	67 nodes to be scored	(559 eliminated genes)
    
    
    	 Level 11:	109 nodes to be scored	(1809 eliminated genes)
    
    
    	 Level 10:	171 nodes to be scored	(3564 eliminated genes)
    
    
    	 Level 9:	243 nodes to be scored	(5410 eliminated genes)
    
    
    	 Level 8:	267 nodes to be scored	(6796 eliminated genes)
    
    
    	 Level 7:	350 nodes to be scored	(9466 eliminated genes)
    
    
    	 Level 6:	362 nodes to be scored	(11754 eliminated genes)
    
    
    	 Level 5:	293 nodes to be scored	(13123 eliminated genes)
    
    
    	 Level 4:	169 nodes to be scored	(14324 eliminated genes)
    
    
    	 Level 3:	82 nodes to be scored	(14897 eliminated genes)
    
    
    	 Level 2:	19 nodes to be scored	(15205 eliminated genes)
    
    
    	 Level 1:	1 nodes to be scored	(15491 eliminated genes)
    



```R
#retrieve the GO enrichment 
goEnrichment   <- GenTable(GOdata, Fisher = results_go, orderBy = "Fisher", topNodes = 100, numChar=1000)
 
```


```R
#graph the GO enrichment 
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
goEnrichment <- goEnrichment[goEnrichment$Fisher < 0.05,] 
goEnrichment <- goEnrichment[goEnrichment$Significant > 2,]#remove any singletons

```


```R
goEnrichment
```


<table class="dataframe">
<caption>A data.frame: 39 × 6</caption>
<thead>
	<tr><th></th><th scope=col>GO.ID</th><th scope=col>Term</th><th scope=col>Annotated</th><th scope=col>Significant</th><th scope=col>Expected</th><th scope=col>Fisher</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0001519</td><td>peptide amidation                                                                </td><td>  19</td><td>12</td><td> 0.22</td><td>2.400e-19</td></tr>
	<tr><th scope=row>2</th><td>GO:0032504</td><td>multicellular organism reproduction                                              </td><td> 969</td><td>24</td><td>11.44</td><td>4.100e-14</td></tr>
	<tr><th scope=row>3</th><td>GO:0009620</td><td>response to fungus                                                               </td><td>  60</td><td>12</td><td> 0.71</td><td>4.100e-13</td></tr>
	<tr><th scope=row>4</th><td>GO:0044719</td><td>regulation of imaginal disc-derived wing size                                    </td><td>  49</td><td>11</td><td> 0.58</td><td>9.000e-12</td></tr>
	<tr><th scope=row>5</th><td>GO:0006508</td><td>proteolysis                                                                      </td><td>1540</td><td>33</td><td>18.17</td><td>1.300e-08</td></tr>
	<tr><th scope=row>6</th><td>GO:0007613</td><td>memory                                                                           </td><td> 109</td><td>11</td><td> 1.29</td><td>6.400e-08</td></tr>
	<tr><th scope=row>7</th><td>GO:0055114</td><td>oxidation-reduction process                                                      </td><td>1003</td><td>34</td><td>11.84</td><td>7.000e-08</td></tr>
	<tr><th scope=row>8</th><td>GO:0030206</td><td>chondroitin sulfate biosynthetic process                                         </td><td>  31</td><td> 5</td><td> 0.37</td><td>2.900e-05</td></tr>
	<tr><th scope=row>10</th><td>GO:2000098</td><td>negative regulation of smooth muscle cell-matrix adhesion                        </td><td>  24</td><td> 4</td><td> 0.28</td><td>1.700e-04</td></tr>
	<tr><th scope=row>11</th><td>GO:0060588</td><td>negative regulation of lipoprotein lipid oxidation                               </td><td>  24</td><td> 4</td><td> 0.28</td><td>1.700e-04</td></tr>
	<tr><th scope=row>12</th><td>GO:2000405</td><td>negative regulation of T cell migration                                          </td><td>  24</td><td> 4</td><td> 0.28</td><td>1.700e-04</td></tr>
	<tr><th scope=row>13</th><td>GO:0018401</td><td>peptidyl-proline hydroxylation to 4-hydroxy-L-proline                            </td><td>  25</td><td> 4</td><td> 0.30</td><td>2.000e-04</td></tr>
	<tr><th scope=row>14</th><td>GO:0071638</td><td>negative regulation of monocyte chemotactic protein-1 production                 </td><td>  27</td><td> 4</td><td> 0.32</td><td>2.700e-04</td></tr>
	<tr><th scope=row>15</th><td>GO:1900016</td><td>negative regulation of cytokine production involved in inflammatory response     </td><td>  28</td><td> 4</td><td> 0.33</td><td>3.100e-04</td></tr>
	<tr><th scope=row>16</th><td>GO:0010642</td><td>negative regulation of platelet-derived growth factor receptor signaling pathway </td><td>  28</td><td> 4</td><td> 0.33</td><td>3.100e-04</td></tr>
	<tr><th scope=row>17</th><td>GO:0048662</td><td>negative regulation of smooth muscle cell proliferation                          </td><td>  31</td><td> 4</td><td> 0.37</td><td>4.600e-04</td></tr>
	<tr><th scope=row>18</th><td>GO:0014012</td><td>peripheral nervous system axon regeneration                                      </td><td>  31</td><td> 4</td><td> 0.37</td><td>4.600e-04</td></tr>
	<tr><th scope=row>19</th><td>GO:0051895</td><td>negative regulation of focal adhesion assembly                                   </td><td>  35</td><td> 4</td><td> 0.41</td><td>7.400e-04</td></tr>
	<tr><th scope=row>20</th><td>GO:0006869</td><td>lipid transport                                                                  </td><td> 346</td><td>13</td><td> 4.08</td><td>8.000e-04</td></tr>
	<tr><th scope=row>21</th><td>GO:0042308</td><td>negative regulation of protein import into nucleus                               </td><td>  43</td><td> 4</td><td> 0.51</td><td>1.620e-03</td></tr>
	<tr><th scope=row>23</th><td>GO:0051592</td><td>response to calcium ion                                                          </td><td>  78</td><td> 5</td><td> 0.92</td><td>2.280e-03</td></tr>
	<tr><th scope=row>24</th><td>GO:0016540</td><td>protein autoprocessing                                                           </td><td>  23</td><td> 3</td><td> 0.27</td><td>2.410e-03</td></tr>
	<tr><th scope=row>29</th><td>GO:0048678</td><td>response to axon injury                                                          </td><td>  70</td><td> 6</td><td> 0.83</td><td>6.860e-03</td></tr>
	<tr><th scope=row>30</th><td>GO:0055085</td><td>transmembrane transport                                                          </td><td>1050</td><td>20</td><td>12.39</td><td>6.870e-03</td></tr>
	<tr><th scope=row>31</th><td>GO:0072347</td><td>response to anesthetic                                                           </td><td>  61</td><td> 5</td><td> 0.72</td><td>6.920e-03</td></tr>
	<tr><th scope=row>32</th><td>GO:0070555</td><td>response to interleukin-1                                                        </td><td>  61</td><td> 3</td><td> 0.72</td><td>7.070e-03</td></tr>
	<tr><th scope=row>33</th><td>GO:0042246</td><td>tissue regeneration                                                              </td><td>  67</td><td> 4</td><td> 0.79</td><td>8.070e-03</td></tr>
	<tr><th scope=row>43</th><td>GO:0016192</td><td>vesicle-mediated transport                                                       </td><td>1340</td><td>16</td><td>15.81</td><td>1.268e-02</td></tr>
	<tr><th scope=row>47</th><td>GO:0007585</td><td>respiratory gaseous exchange                                                     </td><td>  25</td><td> 3</td><td> 0.30</td><td>1.661e-02</td></tr>
	<tr><th scope=row>49</th><td>GO:0030512</td><td>negative regulation of transforming growth factor beta receptor signaling pathway</td><td>  46</td><td> 3</td><td> 0.54</td><td>1.693e-02</td></tr>
	<tr><th scope=row>50</th><td>GO:0008218</td><td>bioluminescence                                                                  </td><td>  46</td><td> 3</td><td> 0.54</td><td>1.693e-02</td></tr>
	<tr><th scope=row>51</th><td>GO:0016051</td><td>carbohydrate biosynthetic process                                                </td><td> 215</td><td> 6</td><td> 2.54</td><td>1.766e-02</td></tr>
	<tr><th scope=row>54</th><td>GO:0009058</td><td>biosynthetic process                                                             </td><td>5743</td><td>32</td><td>67.78</td><td>2.148e-02</td></tr>
	<tr><th scope=row>66</th><td>GO:0006749</td><td>glutathione metabolic process                                                    </td><td> 104</td><td> 5</td><td> 1.23</td><td>2.479e-02</td></tr>
	<tr><th scope=row>68</th><td>GO:0007274</td><td>neuromuscular synaptic transmission                                              </td><td>  58</td><td> 3</td><td> 0.68</td><td>3.107e-02</td></tr>
	<tr><th scope=row>80</th><td>GO:0042391</td><td>regulation of membrane potential                                                 </td><td> 157</td><td> 7</td><td> 1.85</td><td>3.831e-02</td></tr>
	<tr><th scope=row>83</th><td>GO:1903533</td><td>regulation of protein targeting                                                  </td><td> 138</td><td> 5</td><td> 1.63</td><td>4.085e-02</td></tr>
	<tr><th scope=row>85</th><td>GO:0010043</td><td>response to zinc ion                                                             </td><td>  41</td><td> 3</td><td> 0.48</td><td>4.528e-02</td></tr>
	<tr><th scope=row>86</th><td>GO:0009749</td><td>response to glucose                                                              </td><td> 103</td><td> 3</td><td> 1.22</td><td>4.561e-02</td></tr>
</tbody>
</table>




```R
#extract transcript ids that are significantly enriched in the BCN
myterms =goEnrichment$GO.ID 
mygenes = genesInTerm(GOdata, myterms)
```


```R
#extract the transcript ids for each GO term
var=c()
for (i in 1:length(myterms))
{
   myterm <- myterms[i]
   mygenesforterm <- mygenes[myterm][[1]]
   myfactor <- mygenesforterm %in% myInterestingGenes
   mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
   mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
   var[i]=paste(myterm,"genes:",mygenesforterm2)
}

#write.table(var, "BCN_GO_to_mapping.txt", sep="\t", quote=F)
```


```R
#removed GO enrichment terms that have the same exact set of genes 
goEnrichment_subsetted <- read_excel("GO_enrichment_GO_figure_subsetted.xlsx")
```


```R
#plot the GO enrichment results 
ntop = 25
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) 
plot_BCN_GO_BP <- ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = Significant, fill = -log10(Fisher))) +
 ylim(-1,21) + 
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
labs(
    title = 'GO enrichment - BP - BCN')+
  theme_bw(base_size = 20) +
labs(size= "Number of Genes")+
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 13, face = 'bold', vjust = 0.5),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), 
    legend.text = element_text(size = 16, face = "bold"))+


  coord_flip()


```


```R
plot_BCN_GO_BP + labs(x = NULL)
```


![png](output_103_0.png)



```R
#options(repr.plot.width=14, repr.plot.height=8, repr.plot.res = 500)
```

<h2 style="color: #3A9BDC;"> Use topGO to identify enriched molecular functions in BCN  </h2> 


```R
#run the topGO function
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO)

results_go_MF <- runTest(GOdata_MF, algorithm="weight01", statistic="Fisher")


```

    
    Building most specific GOs .....
    
    	( 3583 GO terms found. )
    
    
    Build GO DAG topology ..........
    
    	( 3618 GO terms and 4733 relations. )
    
    
    Annotating nodes ...............
    
    	( 16040 genes annotated to the GO terms. )
    
    
    			 -- Weight01 Algorithm -- 
    
    		 the algorithm is scoring 525 nontrivial nodes
    		 parameters: 
    			 test statistic: fisher
    
    
    	 Level 13:	1 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 12:	8 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 11:	8 nodes to be scored	(35 eliminated genes)
    
    
    	 Level 10:	13 nodes to be scored	(141 eliminated genes)
    
    
    	 Level 9:	24 nodes to be scored	(260 eliminated genes)
    
    
    	 Level 8:	33 nodes to be scored	(405 eliminated genes)
    
    
    	 Level 7:	64 nodes to be scored	(3062 eliminated genes)
    
    
    	 Level 6:	107 nodes to be scored	(3867 eliminated genes)
    
    
    	 Level 5:	104 nodes to be scored	(6554 eliminated genes)
    
    
    	 Level 4:	109 nodes to be scored	(9256 eliminated genes)
    
    
    	 Level 3:	43 nodes to be scored	(13624 eliminated genes)
    
    
    	 Level 2:	10 nodes to be scored	(14595 eliminated genes)
    
    
    	 Level 1:	1 nodes to be scored	(15936 eliminated genes)
    



```R
#retrieve the GO enrichment 
goEnrichment_MF   <- GenTable(GOdata_MF, Fisher = results_go_MF, orderBy = "Fisher", topNodes = 100, numChar=1000)
 
```


```R
goEnrichment_MF$Fisher <- as.numeric(goEnrichment_MF$Fisher)
goEnrichment_MF <- goEnrichment_MF[goEnrichment_MF$Fisher < 0.05,] 
goEnrichment_MF <- goEnrichment_MF[goEnrichment_MF$Significant > 2,]
goEnrichment_MF
```


<table class="dataframe">
<caption>A data.frame: 23 × 6</caption>
<thead>
	<tr><th></th><th scope=col>GO.ID</th><th scope=col>Term</th><th scope=col>Annotated</th><th scope=col>Significant</th><th scope=col>Expected</th><th scope=col>Fisher</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0004598</td><td>peptidylamidoglycolate lyase activity                                                                                        </td><td> 20</td><td>12</td><td>0.24</td><td>1.400e-19</td></tr>
	<tr><th scope=row>2</th><td>GO:0004504</td><td>peptidylglycine monooxygenase activity                                                                                       </td><td> 20</td><td>12</td><td>0.24</td><td>1.400e-19</td></tr>
	<tr><th scope=row>3</th><td>GO:0005507</td><td>copper ion binding                                                                                                           </td><td> 65</td><td>12</td><td>0.79</td><td>2.900e-12</td></tr>
	<tr><th scope=row>4</th><td>GO:0004181</td><td>metallocarboxypeptidase activity                                                                                             </td><td> 36</td><td> 7</td><td>0.44</td><td>7.800e-08</td></tr>
	<tr><th scope=row>5</th><td>GO:0004222</td><td>metalloendopeptidase activity                                                                                                </td><td>118</td><td> 8</td><td>1.44</td><td>3.300e-05</td></tr>
	<tr><th scope=row>6</th><td>GO:0005509</td><td>calcium ion binding                                                                                                          </td><td>640</td><td>19</td><td>7.82</td><td>4.200e-05</td></tr>
	<tr><th scope=row>7</th><td>GO:0046422</td><td>violaxanthin de-epoxidase activity                                                                                           </td><td> 10</td><td> 3</td><td>0.12</td><td>1.300e-04</td></tr>
	<tr><th scope=row>8</th><td>GO:0050659</td><td>N-acetylgalactosamine 4-sulfate 6-O-sulfotransferase activity                                                                </td><td> 12</td><td> 3</td><td>0.15</td><td>2.300e-04</td></tr>
	<tr><th scope=row>9</th><td>GO:0001537</td><td>N-acetylgalactosamine 4-O-sulfotransferase activity                                                                          </td><td> 16</td><td> 3</td><td>0.20</td><td>5.700e-04</td></tr>
	<tr><th scope=row>11</th><td>GO:0004252</td><td>serine-type endopeptidase activity                                                                                           </td><td>241</td><td> 9</td><td>2.94</td><td>9.900e-04</td></tr>
	<tr><th scope=row>13</th><td>GO:0018833</td><td>DDT-dehydrochlorinase activity                                                                                               </td><td> 26</td><td> 3</td><td>0.32</td><td>2.440e-03</td></tr>
	<tr><th scope=row>15</th><td>GO:0004656</td><td>procollagen-proline 4-dioxygenase activity                                                                                   </td><td> 28</td><td> 3</td><td>0.34</td><td>3.020e-03</td></tr>
	<tr><th scope=row>18</th><td>GO:0015485</td><td>cholesterol binding                                                                                                          </td><td> 73</td><td> 4</td><td>0.89</td><td>7.140e-03</td></tr>
	<tr><th scope=row>19</th><td>GO:0031418</td><td>L-ascorbic acid binding                                                                                                      </td><td> 41</td><td> 3</td><td>0.50</td><td>8.910e-03</td></tr>
	<tr><th scope=row>32</th><td>GO:0016757</td><td>transferase activity, transferring glycosyl groups                                                                           </td><td>341</td><td> 5</td><td>4.16</td><td>1.045e-02</td></tr>
	<tr><th scope=row>33</th><td>GO:0047712</td><td>Cypridina-luciferin 2-monooxygenase activity                                                                                 </td><td> 46</td><td> 3</td><td>0.56</td><td>1.222e-02</td></tr>
	<tr><th scope=row>35</th><td>GO:0051287</td><td>NAD binding                                                                                                                  </td><td> 94</td><td> 4</td><td>1.15</td><td>1.692e-02</td></tr>
	<tr><th scope=row>36</th><td>GO:0016702</td><td>oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen</td><td> 63</td><td> 5</td><td>0.77</td><td>1.741e-02</td></tr>
	<tr><th scope=row>37</th><td>GO:0004867</td><td>serine-type endopeptidase inhibitor activity                                                                                 </td><td> 95</td><td> 4</td><td>1.16</td><td>1.752e-02</td></tr>
	<tr><th scope=row>46</th><td>GO:0004364</td><td>glutathione transferase activity                                                                                             </td><td> 61</td><td> 3</td><td>0.74</td><td>2.587e-02</td></tr>
	<tr><th scope=row>58</th><td>GO:0051119</td><td>sugar transmembrane transporter activity                                                                                     </td><td> 66</td><td> 3</td><td>0.81</td><td>3.168e-02</td></tr>
	<tr><th scope=row>59</th><td>GO:0004180</td><td>carboxypeptidase activity                                                                                                    </td><td> 88</td><td>10</td><td>1.07</td><td>3.185e-02</td></tr>
	<tr><th scope=row>60</th><td>GO:0004601</td><td>peroxidase activity                                                                                                          </td><td>105</td><td> 4</td><td>1.28</td><td>4.051e-02</td></tr>
</tbody>
</table>




```R
#plot GO MF enchriment plot 
ntop = 23
ggdata <- goEnrichment_MF[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) 
GO_MF_BCN_plot <- ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = Significant, fill = -log10(Fisher))) +

  ylim(-1,21) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  labs(
    title = 'GO Analysis - MF - BCN')+
   theme_bw(base_size = 20) +
labs(size= "Number of Genes")+
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 13, face = 'bold', vjust = 0.5),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), 
    legend.text = element_text(size = 16, face = "bold"))+


  coord_flip()

```


```R
GO_MF_BCN_plot + labs(x=NULL)
```


![png](output_110_0.png)



```R
#options(repr.plot.width=18, repr.plot.height=8, repr.plot.res = 500)
```

 <h1 style="color: #3A9BDC;"> BCN network connectivity </h1> 
 



WGCNA identifies genes integral to a (module) network using a network characteristic called module membership (Langfelder and Horvath, 2008). Module membership represents connectivity of a gene with other genes within a module and is used to define centralised hub genes (Langfelder and Horvath, 2008). Module membership (MM) is a measure ranging from 0 to 1, with higher values indicating strong connectivity within a module and lower values indicating weak connectivity (Langfelder and Horvath, 2008). Genes with similar expression patterns within a module are not only correlated with the module's overall expression (MM) but also tightly interconnected with other genes in the same module (IC) (Langfelder and Horvath, 2008). In other words, genes that are strongly correlated with the overall expression pattern of their modules (high MM) are also highly connected to other genes within those modules (high IC). There is a  strong correlation between intramodule connectivity and module membership (Langfelder and Horvath, 2008). 



 <h2 style="color: #3A9BDC;">Determine connectivity between intramodular connectivity and module membership   </h2> 
 


```R
#select the same power used for WGCNA analysis in Section 4.3
ADJ1=abs(cor(datExpr,use="p"))^8
Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors)
head(Alldegrees1)

```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>kTotal</th><th scope=col>kWithin</th><th scope=col>kOut</th><th scope=col>kDiff</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_10000_length_2017_cov_3.60449_g7052_i0</th><td>1553.0305073</td><td>1503.6348251</td><td>49.3956822</td><td>1454.23914297</td></tr>
	<tr><th scope=row>NODE_100017_length_140_cov_9.63529_g92797_i0</th><td>  11.1935631</td><td>   3.9747293</td><td> 7.2188338</td><td>  -3.24410443</td></tr>
	<tr><th scope=row>NODE_10001_length_2016_cov_11.7874_g7053_i0</th><td>1409.7578813</td><td>1361.2942003</td><td>48.4636810</td><td>1312.83051935</td></tr>
	<tr><th scope=row>NODE_100020_length_140_cov_8.6_g92800_i0</th><td>   3.6384341</td><td>   0.7672649</td><td> 2.8711692</td><td>  -2.10390425</td></tr>
	<tr><th scope=row>NODE_100025_length_140_cov_5.88235_g92805_i0</th><td>   5.9945145</td><td>   1.4985511</td><td> 4.4959634</td><td>  -2.99741229</td></tr>
	<tr><th scope=row>NODE_100026_length_140_cov_5.69412_g92806_i0</th><td>   0.9903981</td><td>   0.5349826</td><td> 0.4554154</td><td>   0.07956718</td></tr>
</tbody>
</table>




```R
datME=moduleEigengenes(datExpr,dynamicColors)$eigengenes
signif(cor(datME, use="p"), 2)
```


<table class="dataframe">
<caption>A matrix: 14 × 14 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>MEblack</th><th scope=col>MEblue</th><th scope=col>MEbrown</th><th scope=col>MEgreen</th><th scope=col>MEgreenyellow</th><th scope=col>MEgrey</th><th scope=col>MEmagenta</th><th scope=col>MEpink</th><th scope=col>MEpurple</th><th scope=col>MEred</th><th scope=col>MEsalmon</th><th scope=col>MEtan</th><th scope=col>MEturquoise</th><th scope=col>MEyellow</th></tr>
</thead>
<tbody>
	<tr><th scope=row>MEblack</th><td> 1.000</td><td>-0.210</td><td>-0.36</td><td>-0.360</td><td>-0.0430</td><td> 0.1100</td><td>-0.120</td><td>-0.370</td><td>-0.4500</td><td> 0.7000</td><td>-0.058</td><td> 0.210</td><td> 0.057</td><td>-0.17</td></tr>
	<tr><th scope=row>MEblue</th><td>-0.210</td><td> 1.000</td><td>-0.26</td><td>-0.130</td><td>-0.0710</td><td>-0.0890</td><td>-0.069</td><td>-0.310</td><td> 0.4200</td><td>-0.2800</td><td>-0.150</td><td>-0.160</td><td>-0.680</td><td>-0.35</td></tr>
	<tr><th scope=row>MEbrown</th><td>-0.360</td><td>-0.260</td><td> 1.00</td><td> 0.660</td><td> 0.3000</td><td> 0.5200</td><td> 0.410</td><td> 0.450</td><td> 0.4900</td><td>-0.3400</td><td> 0.220</td><td> 0.380</td><td> 0.670</td><td> 0.76</td></tr>
	<tr><th scope=row>MEgreen</th><td>-0.360</td><td>-0.130</td><td> 0.66</td><td> 1.000</td><td> 0.3000</td><td> 0.2900</td><td>-0.045</td><td> 0.300</td><td> 0.2600</td><td>-0.1600</td><td> 0.540</td><td> 0.230</td><td> 0.460</td><td> 0.75</td></tr>
	<tr><th scope=row>MEgreenyellow</th><td>-0.043</td><td>-0.071</td><td> 0.30</td><td> 0.300</td><td> 1.0000</td><td> 0.0130</td><td> 0.130</td><td> 0.150</td><td>-0.0025</td><td> 0.1200</td><td> 0.300</td><td> 0.120</td><td> 0.280</td><td> 0.41</td></tr>
	<tr><th scope=row>MEgrey</th><td> 0.110</td><td>-0.089</td><td> 0.52</td><td> 0.290</td><td> 0.0130</td><td> 1.0000</td><td> 0.150</td><td> 0.420</td><td> 0.3000</td><td>-0.0095</td><td>-0.260</td><td> 0.240</td><td> 0.290</td><td> 0.47</td></tr>
	<tr><th scope=row>MEmagenta</th><td>-0.120</td><td>-0.069</td><td> 0.41</td><td>-0.045</td><td> 0.1300</td><td> 0.1500</td><td> 1.000</td><td> 0.270</td><td> 0.1600</td><td>-0.0810</td><td>-0.190</td><td> 0.027</td><td> 0.300</td><td> 0.20</td></tr>
	<tr><th scope=row>MEpink</th><td>-0.370</td><td>-0.310</td><td> 0.45</td><td> 0.300</td><td> 0.1500</td><td> 0.4200</td><td> 0.270</td><td> 1.000</td><td>-0.0650</td><td>-0.1800</td><td>-0.082</td><td> 0.091</td><td> 0.430</td><td> 0.34</td></tr>
	<tr><th scope=row>MEpurple</th><td>-0.450</td><td> 0.420</td><td> 0.49</td><td> 0.260</td><td>-0.0025</td><td> 0.3000</td><td> 0.160</td><td>-0.065</td><td> 1.0000</td><td>-0.5700</td><td>-0.027</td><td> 0.120</td><td>-0.280</td><td> 0.19</td></tr>
	<tr><th scope=row>MEred</th><td> 0.700</td><td>-0.280</td><td>-0.34</td><td>-0.160</td><td> 0.1200</td><td>-0.0095</td><td>-0.081</td><td>-0.180</td><td>-0.5700</td><td> 1.0000</td><td> 0.110</td><td> 0.160</td><td> 0.170</td><td> 0.12</td></tr>
	<tr><th scope=row>MEsalmon</th><td>-0.058</td><td>-0.150</td><td> 0.22</td><td> 0.540</td><td> 0.3000</td><td>-0.2600</td><td>-0.190</td><td>-0.082</td><td>-0.0270</td><td> 0.1100</td><td> 1.000</td><td> 0.170</td><td> 0.280</td><td> 0.44</td></tr>
	<tr><th scope=row>MEtan</th><td> 0.210</td><td>-0.160</td><td> 0.38</td><td> 0.230</td><td> 0.1200</td><td> 0.2400</td><td> 0.027</td><td> 0.091</td><td> 0.1200</td><td> 0.1600</td><td> 0.170</td><td> 1.000</td><td> 0.340</td><td> 0.27</td></tr>
	<tr><th scope=row>MEturquoise</th><td> 0.057</td><td>-0.680</td><td> 0.67</td><td> 0.460</td><td> 0.2800</td><td> 0.2900</td><td> 0.300</td><td> 0.430</td><td>-0.2800</td><td> 0.1700</td><td> 0.280</td><td> 0.340</td><td> 1.000</td><td> 0.69</td></tr>
	<tr><th scope=row>MEyellow</th><td>-0.170</td><td>-0.350</td><td> 0.76</td><td> 0.750</td><td> 0.4100</td><td> 0.4700</td><td> 0.200</td><td> 0.340</td><td> 0.1900</td><td> 0.1200</td><td> 0.440</td><td> 0.270</td><td> 0.690</td><td> 1.00</td></tr>
</tbody>
</table>




```R
datKME =signedKME(datExpr, datME, outputColumnName="MM.")
head(datKME)
```


<table class="dataframe">
<caption>A data.frame: 6 × 14</caption>
<thead>
	<tr><th></th><th scope=col>MM.black</th><th scope=col>MM.blue</th><th scope=col>MM.brown</th><th scope=col>MM.green</th><th scope=col>MM.greenyellow</th><th scope=col>MM.grey</th><th scope=col>MM.magenta</th><th scope=col>MM.pink</th><th scope=col>MM.purple</th><th scope=col>MM.red</th><th scope=col>MM.salmon</th><th scope=col>MM.tan</th><th scope=col>MM.turquoise</th><th scope=col>MM.yellow</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_10000_length_2017_cov_3.60449_g7052_i0</th><td>-0.15532029</td><td> 0.97601523</td><td>-0.2846444</td><td>-0.144879714</td><td>-0.039132553</td><td>-0.1449490</td><td>-0.0693439986</td><td>-0.34229485</td><td> 0.33625103</td><td>-0.26151146</td><td>-0.11164188</td><td>-0.169304546</td><td>-0.6317741</td><td>-0.379896029</td></tr>
	<tr><th scope=row>NODE_100017_length_140_cov_9.63529_g92797_i0</th><td>-0.26226476</td><td>-0.11715241</td><td> 0.5758819</td><td> 0.571758093</td><td> 0.235266897</td><td> 0.2928844</td><td> 0.2376183428</td><td> 0.31743155</td><td> 0.30416925</td><td>-0.07433324</td><td> 0.28630390</td><td> 0.240246230</td><td> 0.3381483</td><td> 0.572010981</td></tr>
	<tr><th scope=row>NODE_10001_length_2016_cov_11.7874_g7053_i0</th><td>-0.22970154</td><td> 0.95923111</td><td>-0.2430134</td><td>-0.118696037</td><td>-0.051335797</td><td>-0.1268765</td><td>-0.0442088696</td><td>-0.30620975</td><td> 0.37519349</td><td>-0.29075774</td><td>-0.09681686</td><td>-0.111470875</td><td>-0.6219299</td><td>-0.362494605</td></tr>
	<tr><th scope=row>NODE_100020_length_140_cov_8.6_g92800_i0</th><td> 0.23556515</td><td> 0.08400336</td><td>-0.3229051</td><td>-0.283185503</td><td> 0.025327072</td><td>-0.1555593</td><td>-0.0005149578</td><td>-0.44191396</td><td>-0.16374555</td><td> 0.47400625</td><td>-0.11875783</td><td>-0.142910689</td><td>-0.1472292</td><td> 0.002144837</td></tr>
	<tr><th scope=row>NODE_100025_length_140_cov_5.88235_g92805_i0</th><td> 0.52698437</td><td>-0.29706323</td><td> 0.1786959</td><td> 0.062392534</td><td> 0.007372184</td><td> 0.2959813</td><td> 0.0215938414</td><td>-0.08060444</td><td> 0.00623096</td><td> 0.28659794</td><td> 0.11489752</td><td> 0.006724937</td><td> 0.2851930</td><td> 0.300323691</td></tr>
	<tr><th scope=row>NODE_100026_length_140_cov_5.69412_g92806_i0</th><td> 0.02933126</td><td>-0.02126167</td><td> 0.1233537</td><td> 0.008712608</td><td>-0.310629963</td><td> 0.3781600</td><td>-0.0309052488</td><td>-0.17177320</td><td> 0.26500061</td><td>-0.08286882</td><td>-0.16321921</td><td>-0.137215951</td><td>-0.0109511</td><td> 0.077314310</td></tr>
</tbody>
</table>




```R
nrow(datKME)
```


31097



```R
#determine the module membership vs intramodular connectivity for the red module (the BCN)
which.color="red" 
restrictGenes=dynamicColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
(datKME[restrictGenes, paste("MM.", which.color, sep="")])^10,
col=which.color,
xlab="Intramodular Connectivity",
ylab="(Module Membership)^8")
```


![png](output_119_0.png)



```R
#function from the R package repr to visualize figures in Jupyter Notebook(optional)
options(repr.plot.width=8, repr.plot.height=7)
```

 <h2 style="color: #3A9BDC;">Identify genes in the BCN with highest intramodular connectivity </h2> 

Identify all co-expressed transcripts within the BCN (red module) that have a high module membership (MM > 0.8). An MM value of 0.8 or higher indicates a strong association with the module eigengene, suggesting that the gene plays a central role in the module's network and can be considered a hub gene (Langfelder and Horvath, 2008). 



```R
combined_datKME_ADJ1 = merge(datKME, Alldegrees1, by=0)
```


```R
head(combined_datKME_ADJ1)
```


<table class="dataframe">
<caption>A data.frame: 6 × 19</caption>
<thead>
	<tr><th></th><th scope=col>Row.names</th><th scope=col>MM.black</th><th scope=col>MM.blue</th><th scope=col>MM.brown</th><th scope=col>MM.green</th><th scope=col>MM.greenyellow</th><th scope=col>MM.grey</th><th scope=col>MM.magenta</th><th scope=col>MM.pink</th><th scope=col>MM.purple</th><th scope=col>MM.red</th><th scope=col>MM.salmon</th><th scope=col>MM.tan</th><th scope=col>MM.turquoise</th><th scope=col>MM.yellow</th><th scope=col>kTotal</th><th scope=col>kWithin</th><th scope=col>kOut</th><th scope=col>kDiff</th></tr>
	<tr><th></th><th scope=col>&lt;I&lt;chr&gt;&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_1_length_23329_cov_11.1744_g0_i0       </td><td> 0.22666782</td><td> 0.04040074</td><td>-0.19187472</td><td>-0.05986384</td><td>-0.14137351</td><td>-0.15607518</td><td>-0.15311363</td><td>-0.4167876</td><td> 0.05360364</td><td> 0.33473461</td><td> 0.25271613</td><td> 0.1551101</td><td>-0.1389397</td><td> 0.005845331</td><td>   4.438268</td><td>   1.473735</td><td> 2.964533</td><td>  -1.490798</td></tr>
	<tr><th scope=row>2</th><td>NODE_10_length_13445_cov_1.1233_g7_i0       </td><td> 0.19910742</td><td>-0.42991626</td><td> 0.07236738</td><td>-0.03095044</td><td> 0.02067577</td><td> 0.18458297</td><td> 0.17627597</td><td> 0.3534248</td><td>-0.48570092</td><td> 0.31624714</td><td> 0.11926391</td><td> 0.1700352</td><td> 0.4722523</td><td> 0.154050407</td><td>  16.797557</td><td>  11.049569</td><td> 5.747989</td><td>   5.301580</td></tr>
	<tr><th scope=row>3</th><td>NODE_100_length_8069_cov_3.11705_g65_i0     </td><td>-0.09566595</td><td> 0.58874615</td><td>-0.30267560</td><td>-0.11441082</td><td>-0.16633844</td><td>-0.05056483</td><td>-0.13963583</td><td>-0.5034320</td><td> 0.46145430</td><td>-0.14571793</td><td>-0.20310799</td><td>-0.1772821</td><td>-0.7030341</td><td>-0.279742188</td><td> 186.572675</td><td> 154.124679</td><td>32.447996</td><td> 121.676683</td></tr>
	<tr><th scope=row>4</th><td>NODE_10000_length_2017_cov_3.60449_g7052_i0 </td><td>-0.15532029</td><td> 0.97601523</td><td>-0.28464441</td><td>-0.14487971</td><td>-0.03913255</td><td>-0.14494902</td><td>-0.06934400</td><td>-0.3422948</td><td> 0.33625103</td><td>-0.26151146</td><td>-0.11164188</td><td>-0.1693045</td><td>-0.6317741</td><td>-0.379896029</td><td>1553.030507</td><td>1503.634825</td><td>49.395682</td><td>1454.239143</td></tr>
	<tr><th scope=row>5</th><td>NODE_10001_length_2016_cov_11.7874_g7053_i0 </td><td>-0.22970154</td><td> 0.95923111</td><td>-0.24301340</td><td>-0.11869604</td><td>-0.05133580</td><td>-0.12687646</td><td>-0.04420887</td><td>-0.3062098</td><td> 0.37519349</td><td>-0.29075774</td><td>-0.09681686</td><td>-0.1114709</td><td>-0.6219299</td><td>-0.362494605</td><td>1409.757881</td><td>1361.294200</td><td>48.463681</td><td>1312.830519</td></tr>
	<tr><th scope=row>6</th><td>NODE_100017_length_140_cov_9.63529_g92797_i0</td><td>-0.26226476</td><td>-0.11715241</td><td> 0.57588185</td><td> 0.57175809</td><td> 0.23526690</td><td> 0.29288445</td><td> 0.23761834</td><td> 0.3174316</td><td> 0.30416925</td><td>-0.07433324</td><td> 0.28630390</td><td> 0.2402462</td><td> 0.3381483</td><td> 0.572010981</td><td>  11.193563</td><td>   3.974729</td><td> 7.218834</td><td>  -3.244104</td></tr>
</tbody>
</table>




```R
nrow(combined_datKME_ADJ1)
```


31097



```R
subset_combined_datKME_ADJ1 = dplyr::select(combined_datKME_ADJ1, c("Row.names","MM.red", "kWithin","kTotal"))

```


```R
head(subset_combined_datKME_ADJ1)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Row.names</th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th></tr>
	<tr><th></th><th scope=col>&lt;I&lt;chr&gt;&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_1_length_23329_cov_11.1744_g0_i0       </td><td> 0.33473461</td><td>   1.473735</td><td>   4.438268</td></tr>
	<tr><th scope=row>2</th><td>NODE_10_length_13445_cov_1.1233_g7_i0       </td><td> 0.31624714</td><td>  11.049569</td><td>  16.797557</td></tr>
	<tr><th scope=row>3</th><td>NODE_100_length_8069_cov_3.11705_g65_i0     </td><td>-0.14571793</td><td> 154.124679</td><td> 186.572675</td></tr>
	<tr><th scope=row>4</th><td>NODE_10000_length_2017_cov_3.60449_g7052_i0 </td><td>-0.26151146</td><td>1503.634825</td><td>1553.030507</td></tr>
	<tr><th scope=row>5</th><td>NODE_10001_length_2016_cov_11.7874_g7053_i0 </td><td>-0.29075774</td><td>1361.294200</td><td>1409.757881</td></tr>
	<tr><th scope=row>6</th><td>NODE_100017_length_140_cov_9.63529_g92797_i0</td><td>-0.07433324</td><td>   3.974729</td><td>  11.193563</td></tr>
</tbody>
</table>




```R
nrow(subset_combined_datKME_ADJ1)
```


31097



```R
#make the transcript ids rownames again
row.names(subset_combined_datKME_ADJ1) <- subset_combined_datKME_ADJ1$Row.names
subset_combined_datKME_ADJ1[1] <- NULL
```


```R
head(subset_combined_datKME_ADJ1)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_1_length_23329_cov_11.1744_g0_i0</th><td> 0.33473461</td><td>   1.473735</td><td>   4.438268</td></tr>
	<tr><th scope=row>NODE_10_length_13445_cov_1.1233_g7_i0</th><td> 0.31624714</td><td>  11.049569</td><td>  16.797557</td></tr>
	<tr><th scope=row>NODE_100_length_8069_cov_3.11705_g65_i0</th><td>-0.14571793</td><td> 154.124679</td><td> 186.572675</td></tr>
	<tr><th scope=row>NODE_10000_length_2017_cov_3.60449_g7052_i0</th><td>-0.26151146</td><td>1503.634825</td><td>1553.030507</td></tr>
	<tr><th scope=row>NODE_10001_length_2016_cov_11.7874_g7053_i0</th><td>-0.29075774</td><td>1361.294200</td><td>1409.757881</td></tr>
	<tr><th scope=row>NODE_100017_length_140_cov_9.63529_g92797_i0</th><td>-0.07433324</td><td>   3.974729</td><td>  11.193563</td></tr>
</tbody>
</table>




```R
#subset the subset_combined_datKME_ADJ1 to include just transcripts from the red module (the BCN)
top_datKME_ADJ1_.8 <-subset_combined_datKME_ADJ1 %>% dplyr::filter(subset_combined_datKME_ADJ1$MM.red >= 0.79) 

```


```R
head(top_datKME_ADJ1_.8)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_10321_length_1972_cov_1770.41_g3092_i1</th><td>0.8010980</td><td>22.18196</td><td>31.65077</td></tr>
	<tr><th scope=row>NODE_10505_length_1950_cov_293.099_g7393_i0</th><td>0.9121651</td><td>47.49488</td><td>55.68645</td></tr>
	<tr><th scope=row>NODE_10519_length_1949_cov_2.1246_g7404_i0</th><td>0.9043658</td><td>49.06309</td><td>52.70805</td></tr>
	<tr><th scope=row>NODE_106809_length_128_cov_1412.44_g99589_i0</th><td>0.9207360</td><td>81.54745</td><td>89.02335</td></tr>
	<tr><th scope=row>NODE_10748_length_1916_cov_1427.62_g3092_i4</th><td>0.7907177</td><td>20.04188</td><td>28.97377</td></tr>
	<tr><th scope=row>NODE_10872_length_1900_cov_8.46721_g7639_i0</th><td>0.8035119</td><td>20.69150</td><td>26.06230</td></tr>
</tbody>
</table>




```R
#make rowname a column for downstream subsetting 
top_datKME_ADJ1_.8$transcript_id <- rownames(top_datKME_ADJ1_.8)
```


```R
head(top_datKME_ADJ1_.8)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th><th scope=col>transcript_id</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_10321_length_1972_cov_1770.41_g3092_i1</th><td>0.8010980</td><td>22.18196</td><td>31.65077</td><td>NODE_10321_length_1972_cov_1770.41_g3092_i1 </td></tr>
	<tr><th scope=row>NODE_10505_length_1950_cov_293.099_g7393_i0</th><td>0.9121651</td><td>47.49488</td><td>55.68645</td><td>NODE_10505_length_1950_cov_293.099_g7393_i0 </td></tr>
	<tr><th scope=row>NODE_10519_length_1949_cov_2.1246_g7404_i0</th><td>0.9043658</td><td>49.06309</td><td>52.70805</td><td>NODE_10519_length_1949_cov_2.1246_g7404_i0  </td></tr>
	<tr><th scope=row>NODE_106809_length_128_cov_1412.44_g99589_i0</th><td>0.9207360</td><td>81.54745</td><td>89.02335</td><td>NODE_106809_length_128_cov_1412.44_g99589_i0</td></tr>
	<tr><th scope=row>NODE_10748_length_1916_cov_1427.62_g3092_i4</th><td>0.7907177</td><td>20.04188</td><td>28.97377</td><td>NODE_10748_length_1916_cov_1427.62_g3092_i4 </td></tr>
	<tr><th scope=row>NODE_10872_length_1900_cov_8.46721_g7639_i0</th><td>0.8035119</td><td>20.69150</td><td>26.06230</td><td>NODE_10872_length_1900_cov_8.46721_g7639_i0 </td></tr>
</tbody>
</table>




```R
nrow(top_datKME_ADJ1_.8)
```


233



```R
#reorganize the column orders

top_datKME_ADJ1_.8_column_reordered <- top_datKME_ADJ1_.8[, c("transcript_id", "MM.red", "kWithin", "kTotal")] 

head(top_datKME_ADJ1_.8_column_reordered)

```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>transcript_id</th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_10321_length_1972_cov_1770.41_g3092_i1</th><td>NODE_10321_length_1972_cov_1770.41_g3092_i1 </td><td>0.8010980</td><td>22.18196</td><td>31.65077</td></tr>
	<tr><th scope=row>NODE_10505_length_1950_cov_293.099_g7393_i0</th><td>NODE_10505_length_1950_cov_293.099_g7393_i0 </td><td>0.9121651</td><td>47.49488</td><td>55.68645</td></tr>
	<tr><th scope=row>NODE_10519_length_1949_cov_2.1246_g7404_i0</th><td>NODE_10519_length_1949_cov_2.1246_g7404_i0  </td><td>0.9043658</td><td>49.06309</td><td>52.70805</td></tr>
	<tr><th scope=row>NODE_106809_length_128_cov_1412.44_g99589_i0</th><td>NODE_106809_length_128_cov_1412.44_g99589_i0</td><td>0.9207360</td><td>81.54745</td><td>89.02335</td></tr>
	<tr><th scope=row>NODE_10748_length_1916_cov_1427.62_g3092_i4</th><td>NODE_10748_length_1916_cov_1427.62_g3092_i4 </td><td>0.7907177</td><td>20.04188</td><td>28.97377</td></tr>
	<tr><th scope=row>NODE_10872_length_1900_cov_8.46721_g7639_i0</th><td>NODE_10872_length_1900_cov_8.46721_g7639_i0 </td><td>0.8035119</td><td>20.69150</td><td>26.06230</td></tr>
</tbody>
</table>




```R
# reorder spreadsheet based on MM.red column from largest to smallest 

top_datKME_ADJ1_.8_column_reordered_desc <-top_datKME_ADJ1_.8_column_reordered %>% dplyr::arrange(desc(kWithin)) 

```


```R
head(top_datKME_ADJ1_.8_column_reordered_desc)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>transcript_id</th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NODE_37894_length_424_cov_305.016_g30731_i0</th><td>NODE_37894_length_424_cov_305.016_g30731_i0</td><td>0.9395542</td><td>91.07264</td><td>99.28030</td></tr>
	<tr><th scope=row>NODE_35036_length_468_cov_206.31_g27935_i0</th><td>NODE_35036_length_468_cov_206.31_g27935_i0 </td><td>0.9461638</td><td>90.10026</td><td>99.41706</td></tr>
	<tr><th scope=row>NODE_46791_length_332_cov_515.13_g39571_i0</th><td>NODE_46791_length_332_cov_515.13_g39571_i0 </td><td>0.9443668</td><td>89.51240</td><td>99.56021</td></tr>
	<tr><th scope=row>NODE_48716_length_318_cov_339.129_g41496_i0</th><td>NODE_48716_length_318_cov_339.129_g41496_i0</td><td>0.9349035</td><td>88.93315</td><td>95.75253</td></tr>
	<tr><th scope=row>NODE_38913_length_410_cov_355.231_g31740_i0</th><td>NODE_38913_length_410_cov_355.231_g31740_i0</td><td>0.9428738</td><td>87.64126</td><td>95.70440</td></tr>
	<tr><th scope=row>NODE_34295_length_481_cov_68.2207_g27206_i0</th><td>NODE_34295_length_481_cov_68.2207_g27206_i0</td><td>0.9309527</td><td>87.03287</td><td>96.07521</td></tr>
</tbody>
</table>




```R
#annotate the top datKME genes 
top_datKME_ADJ1_.8_column_reordered_desc_annot_trinotate <- left_join(top_datKME_ADJ1_.8_column_reordered_desc,Trinotate_lym_subset,by="transcript_id")

```


```R
head(top_datKME_ADJ1_.8_column_reordered_desc_annot_trinotate)
```


<table class="dataframe">
<caption>A data.frame: 6 × 19</caption>
<thead>
	<tr><th></th><th scope=col>transcript_id</th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>prot_id</th><th scope=col>prot_coords</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_37894_length_424_cov_305.016_g30731_i0</td><td>0.9395542</td><td>91.07264</td><td>99.28030</td><td>g30731</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
	<tr><th scope=row>2</th><td>NODE_35036_length_468_cov_206.31_g27935_i0 </td><td>0.9461638</td><td>90.10026</td><td>99.41706</td><td>g27935</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
	<tr><th scope=row>3</th><td>NODE_46791_length_332_cov_515.13_g39571_i0 </td><td>0.9443668</td><td>89.51240</td><td>99.56021</td><td>g39571</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
	<tr><th scope=row>4</th><td>NODE_48716_length_318_cov_339.129_g41496_i0</td><td>0.9349035</td><td>88.93315</td><td>95.75253</td><td>g41496</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
	<tr><th scope=row>5</th><td>NODE_38913_length_410_cov_355.231_g31740_i0</td><td>0.9428738</td><td>87.64126</td><td>95.70440</td><td>g31740</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
	<tr><th scope=row>6</th><td>NODE_34295_length_481_cov_68.2207_g27206_i0</td><td>0.9309527</td><td>87.03287</td><td>96.07521</td><td>g27206</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#write.csv(top_datKME_ADJ1_.8_column_reordered_desc_annot_trinotate, file = "top_datKME_ADJ1_.8_desc_BCN.csv")
```

<h2 style="color: #3A9BDC;">Identify genes in the bioluminescent upper lip with highest intramodular connectivity </h2> 

Use the results from the DGE analysis in Section Step 9.4 to identify the significantly upregulated genes (expressed uniquely) with the highest module membership (MM > 0.8).


```R
#significantly upregulated genes (expressed uniquely)from section 9.4.4
colnames(unique_genes_bio_upper_lip_info)

```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'transcript_id'</li><li>'baseMean'</li><li>'log2FoldChange'</li><li>'lfcSE'</li><li>'stat'</li><li>'pvalue'</li><li>'padj'</li></ol>




```R
#change the gene column name to transcript ID
colnames(unique_genes_bio_upper_lip_info)[1] <- "transcript_id"
```


```R
#determine how many of the co-expressed genes in the BCN with MM > 0.8 are significantly upregulated in the bioluminescent upper lip

top_datKME_ADJ1_.8_column_reordered_desc_annot_DE_bio_upper_lip <- left_join(top_datKME_ADJ1_.8_column_reordered_desc_annot_trinotate,unique_genes_bio_upper_lip_info,by="transcript_id")

```


```R
head(top_datKME_ADJ1_.8_column_reordered_desc_annot_DE_bio_upper_lip)


```


<table class="dataframe">
<caption>A data.frame: 6 × 25</caption>
<thead>
	<tr><th></th><th scope=col>transcript_id</th><th scope=col>MM.red</th><th scope=col>kWithin</th><th scope=col>kTotal</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>prot_id</th><th scope=col>prot_coords</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>⋯</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_37894_length_424_cov_305.016_g30731_i0</td><td>0.9395542</td><td>91.07264</td><td>99.28030</td><td>g30731</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>2</th><td>NODE_35036_length_468_cov_206.31_g27935_i0 </td><td>0.9461638</td><td>90.10026</td><td>99.41706</td><td>g27935</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>3</th><td>NODE_46791_length_332_cov_515.13_g39571_i0 </td><td>0.9443668</td><td>89.51240</td><td>99.56021</td><td>g39571</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>4</th><td>NODE_48716_length_318_cov_339.129_g41496_i0</td><td>0.9349035</td><td>88.93315</td><td>95.75253</td><td>g41496</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>5</th><td>NODE_38913_length_410_cov_355.231_g31740_i0</td><td>0.9428738</td><td>87.64126</td><td>95.70440</td><td>g31740</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>6</th><td>NODE_34295_length_481_cov_68.2207_g27206_i0</td><td>0.9309527</td><td>87.03287</td><td>96.07521</td><td>g27206</td><td>.</td><td>.</td><td>.</td><td>.</td><td>.</td><td>⋯</td><td>.</td><td>.</td><td>.</td><td>.</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
</tbody>
</table>




```R
#write.csv(top_datKME_ADJ1_.8_column_reordered_desc_annot_DE_bio_upper_lip , file = "top_datKME_ADJ1_.8_annot_DE_unique_Bio_UpperLip.csv")

```

 <h1 style="color: #3A9BDC;"> QC for Differential Gene Expression  </h1> 

Import *V.tsujii* gene expression matrix which consists of three tissue types - upper lip, compound eye and gut - with five biological replicates for each tissue. Generate the sample name sheet (meta sheet) for downstream DESeq2 analyses. 




```R
#read in the gene expression matrix 
organ_level_vt <- read.delim("combined_vargula_tsujii.tab", header = TRUE, sep = "\t", quote = "")
```


```R
head(organ_level_vt)
```


<table class="dataframe">
<caption>A data.frame: 6 × 17</caption>
<thead>
	<tr><th></th><th scope=col>X</th><th scope=col>Vt.1A.counts.tab</th><th scope=col>Vt.1B.counts.tab</th><th scope=col>Vt.1C.counts.tab</th><th scope=col>Vt.2A.counts.tab</th><th scope=col>Vt.2B.counts.tab</th><th scope=col>Vt.2C.counts.tab</th><th scope=col>Vt.3A.counts.tab</th><th scope=col>Vt.3B.counts.tab</th><th scope=col>Vt.3C.counts.tab</th><th scope=col>Vt.4A.counts.tab</th><th scope=col>Vt.4B.counts.tab</th><th scope=col>Vt.4C.counts.tab</th><th scope=col>Vt.5A.counts.tab</th><th scope=col>Vt.5B.counts.tab</th><th scope=col>Vt.5C.counts.tab</th><th scope=col>X.1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;lgl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>NODE_100001_length_141_cov_1.03488_g92781_i0</td><td> 0</td><td>0</td><td> 0</td><td> 0</td><td> 0</td><td>  0</td><td> 0</td><td> 0</td><td>  0</td><td> 0</td><td> 0</td><td> 0</td><td> 0</td><td>0</td><td> 1</td><td>NA</td></tr>
	<tr><th scope=row>2</th><td>NODE_10000_length_2017_cov_3.60449_g7052_i0 </td><td>19</td><td>1</td><td> 8</td><td> 6</td><td> 4</td><td> 81</td><td>40</td><td>15</td><td> 90</td><td>29</td><td> 1</td><td>35</td><td>32</td><td>6</td><td>31</td><td>NA</td></tr>
	<tr><th scope=row>3</th><td>NODE_100017_length_140_cov_9.63529_g92797_i0</td><td> 6</td><td>0</td><td> 2</td><td> 1</td><td> 0</td><td>  5</td><td> 5</td><td> 0</td><td>  9</td><td> 0</td><td> 1</td><td> 7</td><td> 2</td><td>0</td><td> 9</td><td>NA</td></tr>
	<tr><th scope=row>4</th><td>NODE_10001_length_2016_cov_11.7874_g7053_i0 </td><td>23</td><td>2</td><td>29</td><td>48</td><td>20</td><td>180</td><td>66</td><td>23</td><td>220</td><td>40</td><td>30</td><td>51</td><td>49</td><td>2</td><td>85</td><td>NA</td></tr>
	<tr><th scope=row>5</th><td>NODE_100020_length_140_cov_8.6_g92800_i0    </td><td> 0</td><td>0</td><td> 0</td><td> 2</td><td> 0</td><td>  0</td><td> 0</td><td> 0</td><td>  0</td><td> 0</td><td> 0</td><td> 0</td><td> 0</td><td>0</td><td> 0</td><td>NA</td></tr>
	<tr><th scope=row>6</th><td>NODE_100021_length_140_cov_8.04706_g92801_i0</td><td> 0</td><td>0</td><td> 0</td><td> 0</td><td> 0</td><td>  1</td><td> 0</td><td> 0</td><td>  3</td><td> 0</td><td> 0</td><td> 2</td><td> 1</td><td>0</td><td> 0</td><td>NA</td></tr>
</tbody>
</table>




```R
#remove the extra column that got inserted 
row.names(organ_level_vt) <-organ_level_vt$X
organ_level_vt[1]<-NULL

```


```R
organ_level_vt$X.1 <- NULL 
```


```R
#import the sample sheet
sampleinfo <- read_excel("sampleinfo_vtsujii_DGE.xlsx")
```


```R
sampleinfo
```


<table class="dataframe">
<caption>A tibble: 15 × 2</caption>
<thead>
	<tr><th scope=col>samples</th><th scope=col>group</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Vt.1A.counts.tab</td><td>Upper_lip</td></tr>
	<tr><td>Vt.1B.counts.tab</td><td>Eyes     </td></tr>
	<tr><td>Vt.1C.counts.tab</td><td>Gut      </td></tr>
	<tr><td>Vt.2A.counts.tab</td><td>Upper_lip</td></tr>
	<tr><td>Vt.2B.counts.tab</td><td>Eyes     </td></tr>
	<tr><td>Vt.2C.counts.tab</td><td>Gut      </td></tr>
	<tr><td>Vt.3A.counts.tab</td><td>Upper_lip</td></tr>
	<tr><td>Vt.3B.counts.tab</td><td>Eyes     </td></tr>
	<tr><td>Vt.3C.counts.tab</td><td>Gut      </td></tr>
	<tr><td>Vt.4A.counts.tab</td><td>Upper_lip</td></tr>
	<tr><td>Vt.4B.counts.tab</td><td>Eyes     </td></tr>
	<tr><td>Vt.4C.counts.tab</td><td>Gut      </td></tr>
	<tr><td>Vt.5A.counts.tab</td><td>Upper_lip</td></tr>
	<tr><td>Vt.5B.counts.tab</td><td>Eyes     </td></tr>
	<tr><td>Vt.5C.counts.tab</td><td>Gut      </td></tr>
</tbody>
</table>




```R
#generate the DESeq data set 
dds_count_table_organ_level <- DESeqDataSetFromMatrix(countData = organ_level_vt, colData = sampleinfo, design = ~group)

```

    Warning message in DESeqDataSet(se, design = design, ignoreRank):
    “some variables in design formula are characters, converting to factors”



```R
#run DESeq2 function and normalization  
dds_vtsujii <- DESeq(dds_count_table_organ_level, betaPrior = FALSE, parallel = TRUE)

```

    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates: 1 workers
    
    mean-dispersion relationship
    
    final dispersion estimates, fitting model and testing: 1 workers
    



```R
#perform a variance-stabilizing transformation
dds_vtsujii_vsd <- varianceStabilizingTransformation(dds_vtsujii)

```


```R
#transpose the matrix 
sampleDists_dds_vtsujii_vsd <- dist(t(assay(dds_vtsujii_vsd)))

```


```R
#plot the heatmap
sampleDistMatrix_dds_vtsujii_vsd <- as.matrix(sampleDists_dds_vtsujii_vsd)
rownames(sampleDistMatrix_dds_vtsujii_vsd) <- paste(colData(dds_count_table_organ_level)$group) 
colnames(sampleDistMatrix_dds_vtsujii_vsd) <- colData(dds_count_table_organ_level)$samples
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_dds_vtsujii_vsd,
          clustering_distance_rows=sampleDists_dds_vtsujii_vsd,
         clustering_distance_cols=sampleDists_dds_vtsujii_vsd,
         col=colors)

```


![png](output_158_0.png)



```R
options(repr.plot.width=8, repr.plot.height=6, repr.plot.res = 150)
```


```R
pcaData <- plotPCA(dds_vtsujii_vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group, shape=group)) +
labs(color = "Tissue Types")+ labs(shape = "Tissue Types")+
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
   geom_point(size=5) +
theme(panel.grid.major = element_line(colour = "gray97",  size = 1), panel.grid.minor = element_line(linetype = "dotted"), panel.background = element_rect(fill = NA), 
   legend.key = element_rect(fill = "gray100")) + theme(axis.line = element_line(size = 0.5,linetype = "solid")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
coord_fixed() +
scale_color_manual(values = c('#F2C93D','#C97D97','#AC97C9'))+theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 12)  
  )
 
```


![png](output_160_0.png)


 <h1 style="color: #3A9BDC;"> Differential gene expression </h1> 

Differential gene expression analysis was carried out in DESeq2 (Love et al., 2014). For *Vargula tsujii* , we determined differentially upregulated genes in three tissue types - bioluminescent upper lip, compound eye and gut - using five biological replicates for each tissue. DESeq2 was employed using a p-value < 0.05 and FC > 1.5 for the significance of differentially expressed genes using the Benjamini-Hochberg method to account for false discovery rate (FDR). Pairwise comparisons were done across tissue types (i.e., bioluminescent upper lip to compound eye, bioluminescent upper lip to gut, gut to compound eye).  To identify tissue-specific differential expression (i.e. significantly upregulated genes that are uniquely expressed), each tissue was compared to the other two. For example, the expression in the bioluminescent upper lip was determined by comparing it to both the compound eye and the gut tissues.For each pairwise comparison, the reference tissue was specified to determine the significantly upregulated genes in each tissue type (i.e., positive vs negative logfold change). 



```R
#run DESeq2 function and normalization  
dds_vtsujii <- DESeq(dds_count_table_organ_level)

```

    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



<h2 style="color: #3A9BDC;"> DGE - Gut vs Compound Eye  </h2> 



```R
#set Eyes as reference tissue
dds_vtsujii$group <- relevel(dds_vtsujii$group, ref= "Eyes") 
```


```R
#rerun DESeq command after reference is specified
dds_vtsujii <- DESeq(dds_vtsujii)
```

    using pre-existing size factors
    
    estimating dispersions
    
    found already estimated dispersions, replacing these
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



```R
#define contrasts, extract results table, and shrink the log2 fold changes

res_tableOE_unshrunken_Gut_Vs_Eye <- results(dds_vtsujii, contrast= c("group", "Gut", "Eyes") , alpha = 0.05)


res_tableOE_Gut_Vs_Eye  <- lfcShrink(dds_vtsujii, contrast= c("group", "Gut", "Eyes"), res = res_tableOE_unshrunken_Gut_Vs_Eye )


```

    using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
    
    Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
    See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
    Reference: https://doi.org/10.1093/bioinformatics/bty895
    



```R
mcols(res_tableOE_Gut_Vs_Eye, use.names=T)
```


    DataFrame with 6 rows and 2 columns
                           type                               description
                    <character>                               <character>
    baseMean       intermediate mean of normalized counts for all samples
    log2FoldChange      results log2 fold change (MAP): group Gut vs Eyes
    lfcSE               results         standard error: group Gut vs Eyes
    stat                results         Wald statistic: group Gut vs Eyes
    pvalue              results      Wald test p-value: group Gut vs Eyes
    padj                results                      BH adjusted p-values



```R
# set thresholds
# lfc.cutoff value of 0.58 translates to a 1.5 log2 fold change 
# padj.cutoff value of 0.05 

padj.cutoff <- 0.05
lfc.cutoff <- 0.58
```


```R
res_tableOE_Gut_Vs_Eye_tb <- res_tableOE_Gut_Vs_Eye %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
#determine all differentially expressed genes 
sigOE_Gut_Vs_Eye <- res_tableOE_Gut_Vs_Eye_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
```


```R
#extract all genes that are significantly upregulated in the Gut (positive log2 fold change)
sigOE_UPREGULATED_logfold_Gut_vs_Eye <- sigOE_Gut_Vs_Eye %>%
        filter(padj < padj.cutoff & log2FoldChange > lfc.cutoff)

```


```R
head(sigOE_UPREGULATED_logfold_Gut_vs_Eye)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0 </td><td>  7.481691</td><td>4.329005</td><td>1.1986723</td><td>3.646432</td><td>2.659066e-04</td><td>3.236611e-03</td></tr>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>7.152604</td><td>0.9557768</td><td>7.404069</td><td>1.320736e-13</td><td>2.149895e-11</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>5.960708</td><td>1.2577313</td><td>4.683327</td><td>2.822558e-06</td><td>6.652356e-05</td></tr>
	<tr><td>NODE_10047_length_2010_cov_1.26305_g7081_i0 </td><td>  2.108549</td><td>3.092297</td><td>1.1161236</td><td>2.775470</td><td>5.512192e-03</td><td>3.450172e-02</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>6.415144</td><td>1.3379175</td><td>4.675569</td><td>2.931403e-06</td><td>6.869105e-05</td></tr>
	<tr><td>NODE_1006_length_4905_cov_1.0567_g707_i0    </td><td>  2.237450</td><td>3.482954</td><td>1.2248812</td><td>2.814833</td><td>4.880265e-03</td><td>3.151585e-02</td></tr>
</tbody>
</table>




```R
colnames(sigOE_UPREGULATED_logfold_Gut_vs_Eye)[1]<- "transcript_id"
```


```R
#add annotation
sigOE_UPREGULATED_logfold_Gut_vs_Eye_annot <- left_join(sigOE_UPREGULATED_logfold_Gut_vs_Eye,Trinotate_lym_subset,by="transcript_id") 

```


```R
#genes that are significantly upregulated in the Gut (positive log2fold change)
head(sigOE_UPREGULATED_logfold_Gut_vs_Eye_annot)


```


<table class="dataframe">
<caption>A tibble: 6 × 22</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>⋯</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0 </td><td>  7.481691</td><td>4.329005</td><td>1.1986723</td><td>3.646432</td><td>2.659066e-04</td><td>3.236611e-03</td><td>g7057 </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                  </td><td>.                                                                   </td><td>.                  </td><td>.                                                          </td><td>.                                        </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      </td><td>.                                                                                                     </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>7.152604</td><td>0.9557768</td><td>7.404069</td><td>1.320736e-13</td><td>2.149895e-11</td><td>g7060 </td><td>YM9I_CAEEL^YM9I_CAEEL^Q:165-1547,H:47-510^35.789%ID^E:1.6e-90^RecName: Full=Putative serine protease F56F10.1;^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>.</td><td>⋯</td><td>YM9I_CAEEL^YM9I_CAEEL^Q:35-495,H:47-510^35.789%ID^E:7.25e-93^RecName: Full=Putative serine protease F56F10.1;^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                 </td><td>PF05577.12^Peptidase_S28^Serine carboxypeptidase S28^52-486^E:8e-146</td><td>sigP:1^21^0.762^YES</td><td>.                                                          </td><td>ENOG410XSGG^protease, serine, 16 (thymus)</td><td>KEGG:cel:CELE_F56F10.1   </td><td>GO:0045121^cellular_component^membrane raft`GO:0008239^molecular_function^dipeptidyl-peptidase activity`GO:0008236^molecular_function^serine-type peptidase activity`GO:0045087^biological_process^innate immune response`GO:0006508^biological_process^proteolysis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    </td><td>GO:0008236^molecular_function^serine-type peptidase activity`GO:0006508^biological_process^proteolysis</td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>5.960708</td><td>1.2577313</td><td>4.683327</td><td>2.822558e-06</td><td>6.652356e-05</td><td>g7075 </td><td>TSN9_DANRE^TSN9_DANRE^Q:157-810,H:8-230^23.556%ID^E:7.86e-12^RecName: Full=Tetraspanin-9;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii; Teleostei; Ostariophysi; Cypriniformes; Cyprinidae; Danio                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.</td><td>⋯</td><td>TSN9_DANRE^TSN9_DANRE^Q:12-229,H:8-230^26.222%ID^E:3.91e-25^RecName: Full=Tetraspanin-9;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii; Teleostei; Ostariophysi; Cypriniformes; Cyprinidae; Danio                                                                  </td><td>PF00335.20^Tetraspanin^Tetraspanin family^13-226^E:3.8e-43          </td><td>.                  </td><td>ExpAA=90.13^PredHel=4^Topology=i17-39o59-81i88-110o203-225i</td><td>ENOG4111IRY^tetraspanin                  </td><td>KEGG:dre:431733`KO:K17350</td><td>GO:0005887^cellular_component^integral component of plasma membrane                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    </td><td>GO:0016021^cellular_component^integral component of membrane                                          </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10047_length_2010_cov_1.26305_g7081_i0 </td><td>  2.108549</td><td>3.092297</td><td>1.1161236</td><td>2.775470</td><td>5.512192e-03</td><td>3.450172e-02</td><td>g7081 </td><td>DOCK7_HUMAN^DOCK7_HUMAN^Q:1022-1264,H:1900-1980^69.136%ID^E:6.59e-98^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`DOCK7_HUMAN^DOCK7_HUMAN^Q:640-1020,H:1772-1899^53.03%ID^E:6.59e-98^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`DOCK7_HUMAN^DOCK7_HUMAN^Q:1245-1484,H:1974-2053^60%ID^E:6.59e-98^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`DOCK7_HUMAN^DOCK7_HUMAN^Q:1495-1617,H:2056-2096^70.732%ID^E:6.59e-98^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`DOCK7_HUMAN^DOCK7_HUMAN^Q:3-818,H:1561-1830^60.071%ID^E:8.86e-91^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`DOCK7_HUMAN^DOCK7_HUMAN^Q:449-730,H:1712-1801^46.809%ID^E:1.82e-11^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo</td><td>.</td><td>⋯</td><td>DOCK7_HUMAN^DOCK7_HUMAN^Q:1-160,H:1561-1720^83.125%ID^E:1.05e-85^RecName: Full=Dedicator of cytokinesis protein 7;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                      </td><td>PF06920.13^DHR-2^Dock homology region 2^21-164^E:1.1e-46            </td><td>.                  </td><td>.                                                          </td><td>ENOG410XNVY^Dedicator of cytokinesis     </td><td>KEGG:hsa:85440`KO:K21852 </td><td>GO:0030424^cellular_component^axon`GO:0045178^cellular_component^basal part of cell`GO:0005925^cellular_component^focal adhesion`GO:0030426^cellular_component^growth cone`GO:0043005^cellular_component^neuron projection`GO:0005085^molecular_function^guanyl-nucleotide exchange factor activity`GO:0048365^molecular_function^Rac GTPase binding`GO:0090630^biological_process^activation of GTPase activity`GO:0007409^biological_process^axonogenesis`GO:0045200^biological_process^establishment of neuroblast polarity`GO:0022027^biological_process^interkinetic nuclear migration`GO:0000226^biological_process^microtubule cytoskeleton organization`GO:0120163^biological_process^negative regulation of cold-induced thermogenesis`GO:0031175^biological_process^neuron projection development`GO:0033138^biological_process^positive regulation of peptidyl-serine phosphorylation`GO:1904754^biological_process^positive regulation of vascular associated smooth muscle cell migration`GO:0050767^biological_process^regulation of neurogenesis`GO:0007264^biological_process^small GTPase mediated signal transduction</td><td>.                                                                                                     </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>6.415144</td><td>1.3379175</td><td>4.675569</td><td>2.931403e-06</td><td>6.869105e-05</td><td>g93314</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                  </td><td>.                                                                   </td><td>.                  </td><td>.                                                          </td><td>.                                        </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      </td><td>.                                                                                                     </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_1006_length_4905_cov_1.0567_g707_i0    </td><td>  2.237450</td><td>3.482954</td><td>1.2248812</td><td>2.814833</td><td>4.880265e-03</td><td>3.151585e-02</td><td>g707  </td><td>THADA_CANLF^THADA_CANLF^Q:2328-1813,H:1149-1318^33.146%ID^E:4.07e-13^RecName: Full=Thyroid adenoma-associated protein homolog;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Carnivora; Caniformia; Canidae; Canis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>.</td><td>⋯</td><td>THADA_CHLAE^THADA_CHLAE^Q:120-517,H:302-769^19.579%ID^E:2.77e-07^RecName: Full=Thyroid adenoma-associated protein homolog;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Cercopithecidae; Cercopithecinae; Chlorocebus</td><td>.                                                                   </td><td>.                  </td><td>.                                                          </td><td>.                                        </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      </td><td>.                                                                                                     </td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#extract all genes that are significantly upregulated in the Compound Eye (negative log2fold change) but that are downregulated in the Gut. 
sigOE_DOWNREGULATED_logfold_Gut_vs_Eye <- sigOE_Gut_Vs_Eye %>%
        filter(padj < padj.cutoff & log2FoldChange < lfc.cutoff)

```


```R
#save these two dataframes for downstream analysis in Section 9.4

#genes that are significantly upregulated in the Gut (positive log2fold change)
head(sigOE_UPREGULATED_logfold_Gut_vs_Eye)

#genes that are significantly upregulated in the Compound Eye (negative log2fold change) but that are downregulated in the Gut. 
head(sigOE_DOWNREGULATED_logfold_Gut_vs_Eye)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0 </td><td>  7.481691</td><td>4.329005</td><td>1.1986723</td><td>3.646432</td><td>2.659066e-04</td><td>3.236611e-03</td></tr>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>7.152604</td><td>0.9557768</td><td>7.404069</td><td>1.320736e-13</td><td>2.149895e-11</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>5.960708</td><td>1.2577313</td><td>4.683327</td><td>2.822558e-06</td><td>6.652356e-05</td></tr>
	<tr><td>NODE_10047_length_2010_cov_1.26305_g7081_i0 </td><td>  2.108549</td><td>3.092297</td><td>1.1161236</td><td>2.775470</td><td>5.512192e-03</td><td>3.450172e-02</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>6.415144</td><td>1.3379175</td><td>4.675569</td><td>2.931403e-06</td><td>6.869105e-05</td></tr>
	<tr><td>NODE_1006_length_4905_cov_1.0567_g707_i0    </td><td>  2.237450</td><td>3.482954</td><td>1.2248812</td><td>2.814833</td><td>4.880265e-03</td><td>3.151585e-02</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10030_length_2012_cov_88.44_g145_i5    </td><td>77.516326</td><td>-2.004352</td><td>0.7483891</td><td>-2.677266</td><td>7.422565e-03</td><td>4.264394e-02</td></tr>
	<tr><td>NODE_1009_length_4897_cov_2.69228_g709_i0   </td><td>18.842542</td><td>-1.106516</td><td>0.4161180</td><td>-2.659735</td><td>7.820226e-03</td><td>4.417345e-02</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td>32.497969</td><td>-3.786303</td><td>1.2661859</td><td>-3.012530</td><td>2.590799e-03</td><td>1.957288e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td>28.846375</td><td>-4.770626</td><td>0.7722330</td><td>-6.123648</td><td>9.145690e-10</td><td>5.682196e-08</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td> 3.923975</td><td>-3.824371</td><td>1.2262699</td><td>-3.075752</td><td>2.099721e-03</td><td>1.670540e-02</td></tr>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td>15.963167</td><td>-5.219188</td><td>0.9667613</td><td>-5.307557</td><td>1.111044e-07</td><td>3.921827e-06</td></tr>
</tbody>
</table>



<h2 style="color: #3A9BDC;"> DGE - Bioluminescent Upper Lip vs Compound Eye  </h2> 


```R
#define contrasts, extract results table, and shrink the log2 fold changes

res_tableOE_unshrunken_Bio_UpperLip_Vs_Eye <- results(dds_vtsujii, contrast= c("group", "Upper_lip", "Eyes") , alpha = 0.05)


res_tableOE_Bio_UpperLip_Vs_Eye <- lfcShrink(dds_vtsujii, contrast= c("group", "Upper_lip", "Eyes"), res = res_tableOE_unshrunken_Bio_UpperLip_Vs_Eye)



```

    using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
    
    Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
    See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
    Reference: https://doi.org/10.1093/bioinformatics/bty895
    



```R
mcols(res_tableOE_Bio_UpperLip_Vs_Eye, use.names=T)
```


    DataFrame with 6 rows and 2 columns
                           type                                     description
                    <character>                                     <character>
    baseMean       intermediate       mean of normalized counts for all samples
    log2FoldChange      results log2 fold change (MAP): group Upper_lip vs Eyes
    lfcSE               results         standard error: group Upper_lip vs Eyes
    stat                results         Wald statistic: group Upper lip vs Eyes
    pvalue              results      Wald test p-value: group Upper lip vs Eyes
    padj                results                            BH adjusted p-values



```R
res_tableOE_Bio_UpperLip_Vs_Eye_tb <- res_tableOE_Bio_UpperLip_Vs_Eye %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
#determine all differentially expressed genes 
sigOE_Bio_UpperLip_Vs_Eye <- res_tableOE_Bio_UpperLip_Vs_Eye_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
```


```R
#extract all genes that are significantly upregulated in the bioluminescent upper lip (positive log2 fold change)
sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye <- sigOE_Bio_UpperLip_Vs_Eye %>%
        filter(padj < padj.cutoff & log2FoldChange > lfc.cutoff)

```


```R
head(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0</td><td>  7.481691</td><td>4.244836</td><td>1.2023198</td><td>3.570275</td><td>3.566072e-04</td><td>0.0130704556</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1</td><td>144.942119</td><td>3.696833</td><td>0.8140852</td><td>4.529183</td><td>5.921214e-06</td><td>0.0004848792</td></tr>
	<tr><td>NODE_10075_length_2006_cov_50.2255_g7102_i0</td><td> 47.541214</td><td>2.107924</td><td>0.6838391</td><td>3.081562</td><td>2.059176e-03</td><td>0.0439239389</td></tr>
	<tr><td>NODE_10110_length_2001_cov_7.02312_g7128_i0</td><td> 10.296926</td><td>4.431308</td><td>1.0221965</td><td>4.394701</td><td>1.109251e-05</td><td>0.0008282910</td></tr>
	<tr><td>NODE_10314_length_1973_cov_79.1867_g7270_i0</td><td> 51.504733</td><td>1.700903</td><td>0.5199623</td><td>3.269344</td><td>1.077972e-03</td><td>0.0286489771</td></tr>
	<tr><td>NODE_1031_length_4873_cov_1.60834_g724_i0  </td><td> 14.321316</td><td>4.118464</td><td>1.0467162</td><td>3.949132</td><td>7.843501e-05</td><td>0.0040845539</td></tr>
</tbody>
</table>




```R
colnames(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye)[1] <- "transcript_id"
```


```R
#add annotation
sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye_annot <- left_join(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye,Trinotate_lym_subset,by="transcript_id") 

```


```R
#genes that are significantly upregulated in Bioluminescent Upper Lip (positive log2fold change)
head(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye_annot)
```


<table class="dataframe">
<caption>A tibble: 6 × 22</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>⋯</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0</td><td>  7.481691</td><td>4.244836</td><td>1.2023198</td><td>3.570275</td><td>3.566072e-04</td><td>0.0130704556</td><td>g7057</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.</td><td>.</td><td>.                                                                                                                                                                    </td><td>.                       </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>.                                                                                                                                                                                                                                                 </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1</td><td>144.942119</td><td>3.696833</td><td>0.8140852</td><td>4.529183</td><td>5.921214e-06</td><td>0.0004848792</td><td>g4245</td><td>ATS16_MOUSE^ATS16_MOUSE^Q:1751-417,H:93-572^25.662%ID^E:3.7e-36^RecName: Full=A disintegrin and metalloproteinase with thrombospondin motifs 16;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>.</td><td>⋯</td><td>ADT1_CAEEL^ADT1_CAEEL^Q:53-403,H:141-533^27.114%ID^E:4.13e-35^RecName: Full=A disintegrin and metalloproteinase with thrombospondin motifs adt-1 {ECO:0000305};^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 </td><td>PF13688.6^Reprolysin_5^Metallo-peptidase family M12^131-300^E:5.7e-12`PF01421.19^Reprolysin^Reprolysin (M12B) family zinc metalloprotease^135-323^E:1.9e-15`PF13582.6^Reprolysin_3^Metallo-peptidase family M12B Reprolysin-like^142-274^E:4.1e-09`PF13583.6^Reprolysin_4^Metallo-peptidase family M12B Reprolysin-like^200-286^E:3.1e-06`PF13574.6^Reprolysin_2^Metallo-peptidase family M12B Reprolysin-like^219-311^E:3.9e-08`PF17771.1^ADAM_CR_2^ADAM cysteine-rich domain^339-403^E:1.2e-09`PF17771.1^ADAM_CR_2^ADAM cysteine-rich domain^430-498^E:5.6e-05</td><td>.</td><td>.</td><td>ENOG41104P0^Thrombospondin type 1 domain                                                                                                                             </td><td>KEGG:cel:CELE_C02B4.1   </td><td>GO:0005576^cellular_component^extracellular region`GO:0046872^molecular_function^metal ion binding`GO:0004222^molecular_function^metalloendopeptidase activity                                                                                                                                                                                                                                                                                                                                                                                                                   </td><td>GO:0004222^molecular_function^metalloendopeptidase activity`GO:0006508^biological_process^proteolysis                                                                                                                                             </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10075_length_2006_cov_50.2255_g7102_i0</td><td> 47.541214</td><td>2.107924</td><td>0.6838391</td><td>3.081562</td><td>2.059176e-03</td><td>0.0439239389</td><td>g7102</td><td>6PGD_HUMAN^6PGD_HUMAN^Q:1639-203,H:3-483^76.923%ID^E:0^RecName: Full=6-phosphogluconate dehydrogenase, decarboxylating;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      </td><td>.</td><td>⋯</td><td>6PGD_HUMAN^6PGD_HUMAN^Q:1-468,H:14-483^76.596%ID^E:0^RecName: Full=6-phosphogluconate dehydrogenase, decarboxylating;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>PF03446.15^NAD_binding_2^NAD binding domain of 6-phosphogluconate dehydrogenase^1-158^E:3.6e-45`PF00393.19^6PGD^6-phosphogluconate dehydrogenase, C-terminal domain^166-454^E:3.7e-132                                                                                                                                                                                                                                                                                                                                                                          </td><td>.</td><td>.</td><td>COG0362^Catalyzes the oxidative decarboxylation of 6- phosphogluconate to ribulose 5-phosphate and CO(2), with concomitant reduction of NADP to NADPH (By similarity)</td><td>KEGG:hsa:5226`KO:K00033 </td><td>GO:0005829^cellular_component^cytosol`GO:0070062^cellular_component^extracellular exosome`GO:0005634^cellular_component^nucleus`GO:0050661^molecular_function^NADP binding`GO:0004616^molecular_function^phosphogluconate dehydrogenase (decarboxylating) activity`GO:0046177^biological_process^D-gluconate catabolic process`GO:0055114^biological_process^oxidation-reduction process`GO:0019322^biological_process^pentose biosynthetic process`GO:0006098^biological_process^pentose-phosphate shunt`GO:0009051^biological_process^pentose-phosphate shunt, oxidative branch</td><td>GO:0050661^molecular_function^NADP binding`GO:0004616^molecular_function^phosphogluconate dehydrogenase (decarboxylating) activity`GO:0006098^biological_process^pentose-phosphate shunt`GO:0055114^biological_process^oxidation-reduction process</td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10110_length_2001_cov_7.02312_g7128_i0</td><td> 10.296926</td><td>4.431308</td><td>1.0221965</td><td>4.394701</td><td>1.109251e-05</td><td>0.0008282910</td><td>g7128</td><td>GORAB_BOVIN^GORAB_BOVIN^Q:1899-961,H:1-293^32.71%ID^E:5.24e-28^RecName: Full=RAB6-interacting golgin;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  </td><td>.</td><td>⋯</td><td>GORAB_BOVIN^GORAB_BOVIN^Q:1-313,H:1-293^33.028%ID^E:2.23e-39^RecName: Full=RAB6-interacting golgin;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Laurasiatheria; Cetartiodactyla; Ruminantia; Pecora; Bovidae; Bovinae; Bos                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 </td><td>PF04949.13^Transcrip_act^Transcriptional activator^152-296^E:2.8e-06                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>.</td><td>.</td><td>ENOG4111VTM^Golgin, RAB6-interacting                                                                                                                                 </td><td>.                       </td><td>GO:0005794^cellular_component^Golgi apparatus`GO:1905515^biological_process^non-motile cilium assembly                                                                                                                                                                                                                                                                                                                                                                                                                                                                           </td><td>.                                                                                                                                                                                                                                                 </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10314_length_1973_cov_79.1867_g7270_i0</td><td> 51.504733</td><td>1.700903</td><td>0.5199623</td><td>3.269344</td><td>1.077972e-03</td><td>0.0286489771</td><td>g7270</td><td>NSA2_HUMAN^NSA2_HUMAN^Q:1874-1107,H:1-256^77.344%ID^E:1.38e-136^RecName: Full=Ribosome biogenesis protein NSA2 homolog;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      </td><td>.</td><td>⋯</td><td>NSA2_HUMAN^NSA2_HUMAN^Q:1-260,H:1-260^77.692%ID^E:2.05e-153^RecName: Full=Ribosome biogenesis protein NSA2 homolog;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       </td><td>PF01201.22^Ribosomal_S8e^Ribosomal protein S8e^36-259^E:5e-48                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   </td><td>.</td><td>.</td><td>COG2007^40S ribosomal protein S8                                                                                                                                     </td><td>KEGG:hsa:10412`KO:K14842</td><td>GO:0005730^cellular_component^nucleolus`GO:0030687^cellular_component^preribosome, large subunit precursor`GO:0003723^molecular_function^RNA binding`GO:0006364^biological_process^rRNA processing                                                                                                                                                                                                                                                                                                                                                                               </td><td>.                                                                                                                                                                                                                                                 </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_1031_length_4873_cov_1.60834_g724_i0  </td><td> 14.321316</td><td>4.118464</td><td>1.0467162</td><td>3.949132</td><td>7.843501e-05</td><td>0.0040845539</td><td>g724 </td><td>ZN710_HUMAN^ZN710_HUMAN^Q:1685-2203,H:262-427^32.768%ID^E:8.13e-16^RecName: Full=Zinc finger protein 710;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`ZN710_HUMAN^ZN710_HUMAN^Q:1895-2200,H:493-594^31.068%ID^E:2.92e-08^RecName: Full=Zinc finger protein 710;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`ZN710_HUMAN^ZN710_HUMAN^Q:1895-2197,H:381-481^33.333%ID^E:9.07e-08^RecName: Full=Zinc finger protein 710;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`ZN710_HUMAN^ZN710_HUMAN^Q:1895-2200,H:465-566^29.126%ID^E:1.62e-07^RecName: Full=Zinc finger protein 710;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo`ZN710_HUMAN^ZN710_HUMAN^Q:1829-2200,H:389-510^28.906%ID^E:2.63e-07^RecName: Full=Zinc finger protein 710;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo</td><td>.</td><td>⋯</td><td>ZG62_XENLA^ZG62_XENLA^Q:548-653,H:8-113^38.318%ID^E:2.47e-17^RecName: Full=Gastrula zinc finger protein XlCGF62.1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Amphibia; Batrachia; Anura; Pipoidea; Pipidae; Xenopodinae; Xenopus; Xenopus`ZG62_XENLA^ZG62_XENLA^Q:548-642,H:64-159^35.417%ID^E:1.84e-12^RecName: Full=Gastrula zinc finger protein XlCGF62.1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Amphibia; Batrachia; Anura; Pipoidea; Pipidae; Xenopodinae; Xenopus; Xenopus`ZG62_XENLA^ZG62_XENLA^Q:571-648,H:4-80^38.462%ID^E:2.62e-10^RecName: Full=Gastrula zinc finger protein XlCGF62.1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Amphibia; Batrachia; Anura; Pipoidea; Pipidae; Xenopodinae; Xenopus; Xenopus</td><td>PF00096.26^zf-C2H2^Zinc finger, C2H2 type^573-596^E:0.00011`PF00096.26^zf-C2H2^Zinc finger, C2H2 type^602-624^E:0.00036                                                                                                                                                                                                                                                                                                                                                                                                                                         </td><td>.</td><td>.</td><td>.                                                                                                                                                                    </td><td>.                       </td><td>GO:0005634^cellular_component^nucleus`GO:0003677^molecular_function^DNA binding`GO:0046872^molecular_function^metal ion binding                                                                                                                                                                                                                                                                                                                                                                                                                                                  </td><td>GO:0003676^molecular_function^nucleic acid binding                                                                                                                                                                                                </td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#extract all genes that are significantly upregulated in the Compound Eye (negative log2fold change) but that are downregulated in the bioluminescent upper lip 
sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye <- sigOE_Bio_UpperLip_Vs_Eye %>%
        filter(padj < padj.cutoff & log2FoldChange < lfc.cutoff)

```


```R
colnames(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye)[1] <- "transcript_id"
```


```R
#add annotation 
sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye_annot <- left_join(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye,Trinotate_lym_subset,by="transcript_id") 

```


```R
head(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye_annot)
```


<table class="dataframe">
<caption>A tibble: 6 × 22</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>⋯</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td> 15.963167</td><td>-3.442293</td><td>0.9226241</td><td> -3.704888</td><td>2.114840e-04</td><td>8.941592e-03</td><td>g7138 </td><td>GST1D_ANOGA^GST1D_ANOGA^Q:1957-1379,H:1-193^50.515%ID^E:1.3e-57^RecName: Full=Glutathione S-transferase 1, isoform D;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Holometabola; Diptera; Nematocera; Culicoidea; Culicidae; Anophelinae; Anopheles</td><td>.</td><td>⋯</td><td>GST1D_ANOGA^GST1D_ANOGA^Q:3-191,H:1-189^51.579%ID^E:1.46e-69^RecName: Full=Glutathione S-transferase 1, isoform D;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Holometabola; Diptera; Nematocera; Culicoidea; Culicidae; Anophelinae; Anopheles</td><td>PF13417.6^GST_N_3^Glutathione S-transferase, N-terminal domain^5-79^E:6.6e-13`PF02798.20^GST_N^Glutathione S-transferase, N-terminal domain^15-74^E:4.2e-11`PF13409.6^GST_N_2^Glutathione S-transferase, N-terminal domain^37-76^E:3.2e-06`PF00043.25^GST_C^Glutathione S-transferase, C-terminal domain^120-188^E:1.6e-05`PF14497.6^GST_C_3^Glutathione S-transferase, C-terminal domain^123-191^E:6e-06</td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>GO:0005737^cellular_component^cytoplasm`GO:0018833^molecular_function^DDT-dehydrochlorinase activity`GO:0004364^molecular_function^glutathione transferase activity`GO:0006749^biological_process^glutathione metabolic process                                                                                                                                       </td><td>GO:0005515^molecular_function^protein binding                          </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_101633_length_138_cov_2.06024_g94413_i0</td><td>  7.667808</td><td>-6.149554</td><td>1.0781402</td><td> -5.516364</td><td>3.460855e-08</td><td>4.918433e-06</td><td>g94413</td><td>.                                                                                                                                                                                                                                                                                        </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                        </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10218_length_1985_cov_18.5554_g7198_i0 </td><td> 61.993605</td><td>-6.085923</td><td>0.6771456</td><td> -8.841580</td><td>9.437437e-19</td><td>3.922413e-16</td><td>g7198 </td><td>WSD_ACIAD^WSD_ACIAD^Q:1368-265,H:62-446^21.76%ID^E:1.15e-12^RecName: Full=O-acyltransferase WSD;^Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter                                                                                            </td><td>.</td><td>⋯</td><td>WSD_ACIAD^WSD_ACIAD^Q:196-563,H:62-446^22.005%ID^E:9.38e-13^RecName: Full=O-acyltransferase WSD;^Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter                                                                                         </td><td>PF06974.13^DUF1298^Protein of unknown function (DUF1298)^418-558^E:2.4e-27                                                                                                                                                                                                                                                                                                                               </td><td>.                  </td><td>ExpAA=83.94^PredHel=4^Topology=i36-58o73-95i452-474o504-526i</td><td>ENOG410XS7A^Acyltransferase, ws dgat mgat</td><td>KEGG:aci:ACIAD0832`KO:K00635</td><td>GO:0102966^molecular_function^arachidoyl-CoA:1-dodecanol O-acyltransferase activity`GO:0004144^molecular_function^diacylglycerol O-acyltransferase activity`GO:0047196^molecular_function^long-chain-alcohol O-fatty-acyltransferase activity`GO:0006071^biological_process^glycerol metabolic process`GO:0019432^biological_process^triglyceride biosynthetic process</td><td>GO:0004144^molecular_function^diacylglycerol O-acyltransferase activity</td><td>.</td><td>.</td></tr>
	<tr><td>NODE_1021_length_4884_cov_10.72_g718_i0     </td><td> 45.575618</td><td>-6.514480</td><td>0.8450182</td><td> -7.517323</td><td>5.590900e-14</td><td>1.734596e-11</td><td>g718  </td><td>.                                                                                                                                                                                                                                                                                        </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                        </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10232_length_1983_cov_42.0814_g7210_i0 </td><td>434.539937</td><td>-9.286811</td><td>0.7189709</td><td>-12.613262</td><td>1.784501e-36</td><td>6.551496e-33</td><td>g7210 </td><td>.                                                                                                                                                                                                                                                                                        </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                        </td><td>sigP:1^19^0.705^YES</td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_102471_length_136_cov_2.41975_g95251_i0</td><td>  3.353358</td><td>-4.786748</td><td>1.3410759</td><td> -3.571766</td><td>3.545819e-04</td><td>1.301788e-02</td><td>g95251</td><td>.                                                                                                                                                                                                                                                                                        </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                        </td><td>.                  </td><td>.                                                           </td><td>.                                        </td><td>.                           </td><td>.                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                      </td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#save these two dataframes for downstream analysis in Section 9.4

#genes that are significantly upregulated in Bioluminescent Upper Lip (positive log2fold change) 
head(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye)

#genes that are significantly upregulated in the Compound Eye (negative log2fold change) but that are downregulated in the bioluminescent upper lip 
head(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0</td><td>  7.481691</td><td>4.244836</td><td>1.2023198</td><td>3.570275</td><td>3.566072e-04</td><td>0.0130704556</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1</td><td>144.942119</td><td>3.696833</td><td>0.8140852</td><td>4.529183</td><td>5.921214e-06</td><td>0.0004848792</td></tr>
	<tr><td>NODE_10075_length_2006_cov_50.2255_g7102_i0</td><td> 47.541214</td><td>2.107924</td><td>0.6838391</td><td>3.081562</td><td>2.059176e-03</td><td>0.0439239389</td></tr>
	<tr><td>NODE_10110_length_2001_cov_7.02312_g7128_i0</td><td> 10.296926</td><td>4.431308</td><td>1.0221965</td><td>4.394701</td><td>1.109251e-05</td><td>0.0008282910</td></tr>
	<tr><td>NODE_10314_length_1973_cov_79.1867_g7270_i0</td><td> 51.504733</td><td>1.700903</td><td>0.5199623</td><td>3.269344</td><td>1.077972e-03</td><td>0.0286489771</td></tr>
	<tr><td>NODE_1031_length_4873_cov_1.60834_g724_i0  </td><td> 14.321316</td><td>4.118464</td><td>1.0467162</td><td>3.949132</td><td>7.843501e-05</td><td>0.0040845539</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td> 15.963167</td><td>-3.442293</td><td>0.9226241</td><td> -3.704888</td><td>2.114840e-04</td><td>8.941592e-03</td></tr>
	<tr><td>NODE_101633_length_138_cov_2.06024_g94413_i0</td><td>  7.667808</td><td>-6.149554</td><td>1.0781402</td><td> -5.516364</td><td>3.460855e-08</td><td>4.918433e-06</td></tr>
	<tr><td>NODE_10218_length_1985_cov_18.5554_g7198_i0 </td><td> 61.993605</td><td>-6.085923</td><td>0.6771456</td><td> -8.841580</td><td>9.437437e-19</td><td>3.922413e-16</td></tr>
	<tr><td>NODE_1021_length_4884_cov_10.72_g718_i0     </td><td> 45.575618</td><td>-6.514480</td><td>0.8450182</td><td> -7.517323</td><td>5.590900e-14</td><td>1.734596e-11</td></tr>
	<tr><td>NODE_10232_length_1983_cov_42.0814_g7210_i0 </td><td>434.539937</td><td>-9.286811</td><td>0.7189709</td><td>-12.613262</td><td>1.784501e-36</td><td>6.551496e-33</td></tr>
	<tr><td>NODE_102471_length_136_cov_2.41975_g95251_i0</td><td>  3.353358</td><td>-4.786748</td><td>1.3410759</td><td> -3.571766</td><td>3.545819e-04</td><td>1.301788e-02</td></tr>
</tbody>
</table>



 <h2 style="color: #3A9BDC;"> DGE - Bioluminescent Upper Lip vs Gut </h2> 


```R
#now set gut as reference
dds_vtsujii$group <- relevel(dds_vtsujii$group, ref= "Gut") 
```


```R
#rerun DESeq after setting a new reference
dds_vtsujii <- DESeq(dds_vtsujii)
```

    using pre-existing size factors
    
    estimating dispersions
    
    found already estimated dispersions, replacing these
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



```R
#Define contrasts, extract results table, and shrink the log2 fold changes

res_tableOE_unshrunken_Bio_Upper_Lip_Vs_Gut <- results(dds_vtsujii, contrast= c("group", "Upper_lip", "Gut") , alpha = 0.05)


res_tableOE_Bio_Upper_Lip_Vs_Gut  <- lfcShrink(dds_vtsujii, contrast= c("group", "Upper_lip", "Gut"), res = res_tableOE_unshrunken_Bio_Upper_Lip_Vs_Gut )


```

    using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
    
    Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
    See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
    Reference: https://doi.org/10.1093/bioinformatics/bty895
    



```R
mcols(res_tableOE_Bio_Upper_Lip_Vs_Gut, use.names=T)
```


    DataFrame with 6 rows and 2 columns
                           type                                    description
                    <character>                                    <character>
    baseMean       intermediate      mean of normalized counts for all samples
    log2FoldChange      results log2 fold change (MAP): group Upper_lip vs Gut
    lfcSE               results         standard error: group Upper_lip vs Gut
    stat                results         Wald statistic: group Upper lip vs Gut
    pvalue              results      Wald test p-value: group Upper lip vs Gut
    padj                results                           BH adjusted p-values



```R
res_tableOE_Bio_Upper_Lip_Vs_Gut_tb <- res_tableOE_Bio_Upper_Lip_Vs_Gut %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
#all differentially expressed genes 
sigOE_Bio_Upper_Lip_Vs_Gut <- res_tableOE_Bio_Upper_Lip_Vs_Gut_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
```


```R
#extract all genes that are significantly upregulated in the bioluminescent upper lip (positive log2 fold change)
sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut <- sigOE_Bio_Upper_Lip_Vs_Gut %>%
        filter(padj < padj.cutoff & log2FoldChange > lfc.cutoff)

```


```R
head(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td></tr>
</tbody>
</table>




```R
colnames(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)[1] <- "transcript_id"

```


```R
#add annotation 
sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut_annot <- left_join(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut,Trinotate_lym_subset,by="transcript_id") 

```


```R
head(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut_annot)
```


<table class="dataframe">
<caption>A tibble: 6 × 22</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>⋯</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td><td>g7072 </td><td>.                                                                                                                                                                                                                                                                                                                         </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.                  </td><td>.</td><td>.                                       </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                    </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td><td>g4245 </td><td>ATS16_MOUSE^ATS16_MOUSE^Q:1751-417,H:93-572^25.662%ID^E:3.7e-36^RecName: Full=A disintegrin and metalloproteinase with thrombospondin motifs 16;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus</td><td>.</td><td>⋯</td><td>ADT1_CAEEL^ADT1_CAEEL^Q:53-403,H:141-533^27.114%ID^E:4.13e-35^RecName: Full=A disintegrin and metalloproteinase with thrombospondin motifs adt-1 {ECO:0000305};^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          </td><td>PF13688.6^Reprolysin_5^Metallo-peptidase family M12^131-300^E:5.7e-12`PF01421.19^Reprolysin^Reprolysin (M12B) family zinc metalloprotease^135-323^E:1.9e-15`PF13582.6^Reprolysin_3^Metallo-peptidase family M12B Reprolysin-like^142-274^E:4.1e-09`PF13583.6^Reprolysin_4^Metallo-peptidase family M12B Reprolysin-like^200-286^E:3.1e-06`PF13574.6^Reprolysin_2^Metallo-peptidase family M12B Reprolysin-like^219-311^E:3.9e-08`PF17771.1^ADAM_CR_2^ADAM cysteine-rich domain^339-403^E:1.2e-09`PF17771.1^ADAM_CR_2^ADAM cysteine-rich domain^430-498^E:5.6e-05</td><td>.                  </td><td>.</td><td>ENOG41104P0^Thrombospondin type 1 domain</td><td>KEGG:cel:CELE_C02B4.1    </td><td>GO:0005576^cellular_component^extracellular region`GO:0046872^molecular_function^metal ion binding`GO:0004222^molecular_function^metalloendopeptidase activity                                                                                                                                                                                                                                                                                                                        </td><td>GO:0004222^molecular_function^metalloendopeptidase activity`GO:0006508^biological_process^proteolysis</td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td><td>g7115 </td><td>.                                                                                                                                                                                                                                                                                                                         </td><td>.</td><td>⋯</td><td>LRP1_CHICK^LRP1_CHICK^Q:140-221,H:978-1057^56.098%ID^E:4.41e-19^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-220,H:852-936^47.059%ID^E:1.71e-16^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:3572-3641^52.778%ID^E:7.37e-15^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:3491-3564^55.405%ID^E:1.68e-14^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:118-220,H:3514-3624^41.071%ID^E:3.52e-14^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:934-1006^53.425%ID^E:1.45e-13^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-224,H:3610-3697^44.944%ID^E:2.73e-12^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-209,H:893-966^48.649%ID^E:3.5e-12^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:136-208,H:3330-3402^49.315%ID^E:5.14e-12^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-220,H:3450-3535^45.349%ID^E:1.51e-11^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-212,H:2732-2808^51.282%ID^E:1.73e-11^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:1013-1092^48.75%ID^E:3.75e-11^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-209,H:2690-2764^45.57%ID^E:1.5e-10^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:128-220,H:3363-3452^40.206%ID^E:2.28e-09^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-238,H:29-146^36.975%ID^E:1.22e-08^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:2856-2933^45%ID^E:2.14e-08^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:3410-3482^50%ID^E:2.73e-08^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:142-242,H:3698-3804^38.532%ID^E:4.25e-07^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:141-245,H:1107-1222^35.537%ID^E:1.24e-06^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus`LRP1_CHICK^LRP1_CHICK^Q:137-208,H:2560-2630^42.466%ID^E:2.15e-06^RecName: Full=Low-density lipoprotein receptor-related protein 1;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Galloanserae; Galliformes; Phasianidae; Phasianinae; Gallus</td><td>PF00057.18^Ldl_recept_a^Low-density lipoprotein receptor domain class A^136-170^E:2e-09`PF00057.18^Ldl_recept_a^Low-density lipoprotein receptor domain class A^176-209^E:2.5e-11`PF00059.21^Lectin_C^Lectin C-type domain^486-595^E:5.4e-13                                                                                                                                                                                                                                                                                                                    </td><td>sigP:1^23^0.647^YES</td><td>.</td><td>.                                       </td><td>KEGG:gga:396170`KO:K04550</td><td>GO:0005905^cellular_component^clathrin-coated pit`GO:0016021^cellular_component^integral component of membrane`GO:0016964^molecular_function^alpha-2 macroglobulin receptor activity`GO:0005509^molecular_function^calcium ion binding                                                                                                                                                                                                                                                </td><td>GO:0005515^molecular_function^protein binding                                                        </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td><td>g7125 </td><td>.                                                                                                                                                                                                                                                                                                                         </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.                  </td><td>.</td><td>.                                       </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                    </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td><td>g7126 </td><td>MRC2_HUMAN^MRC2_HUMAN^Q:318-1031,H:684-948^22.101%ID^E:1.19e-08^RecName: Full=C-type mannose receptor 2;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                       </td><td>.</td><td>⋯</td><td>MRC2_HUMAN^MRC2_HUMAN^Q:49-286,H:684-948^22.101%ID^E:9.87e-10^RecName: Full=C-type mannose receptor 2;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>PF00059.21^Lectin_C^Lectin C-type domain^62-167^E:1.4e-14`PF00059.21^Lectin_C^Lectin C-type domain^184-287^E:2.8e-08                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>sigP:1^23^0.485^YES</td><td>.</td><td>ENOG410XQ89^Mannose receptor, C type 2  </td><td>KEGG:hsa:9902`KO:K06560  </td><td>GO:0005925^cellular_component^focal adhesion`GO:0016021^cellular_component^integral component of membrane`GO:0016020^cellular_component^membrane`GO:0030246^molecular_function^carbohydrate binding`GO:0005518^molecular_function^collagen binding`GO:0004888^molecular_function^transmembrane signaling receptor activity`GO:0030574^biological_process^collagen catabolic process`GO:0006897^biological_process^endocytosis`GO:0001649^biological_process^osteoblast differentiation</td><td>.                                                                                                    </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td><td>g93955</td><td>.                                                                                                                                                                                                                                                                                                                         </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.                  </td><td>.</td><td>.                                       </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>.                                                                                                    </td><td>.</td><td>.</td></tr>
</tbody>
</table>




```R
#extract all genes that are significantly upregulated in the Gut (negative log2fold change) but that are downregulated in the bioluminescent upper lip. 
sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut <- sigOE_Bio_Upper_Lip_Vs_Gut %>%
        filter(padj < padj.cutoff & log2FoldChange < lfc.cutoff)
```


```R
colnames(sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)[1]<- "transcript_id"
```


```R
#add annotation
sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut_annot <- left_join(sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut, Trinotate_lym_subset,by="transcript_id")
```


```R
#save these two dataframes for downstream analysis in Section 9.3

#genes that are significantly upregulated in the bioluminescent upper lip (positive log2 fold change)
head(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)

#genes that are significantly upregulated in the gut (negative log2fold change) but that are downregulated in the bioluminescent upper lip. 
head(sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut_annot)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 22</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>#gene_id</th><th scope=col>sprot_Top_BLASTX_hit</th><th scope=col>RNAMMER</th><th scope=col>⋯</th><th scope=col>sprot_Top_BLASTP_hit</th><th scope=col>Pfam</th><th scope=col>SignalP</th><th scope=col>TmHMM</th><th scope=col>eggnog</th><th scope=col>Kegg</th><th scope=col>gene_ontology_blast</th><th scope=col>gene_ontology_pfam</th><th scope=col>transcript</th><th scope=col>peptide</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>-6.417432</td><td>0.9300257</td><td>-6.864408</td><td>6.676754e-12</td><td>4.486477e-10</td><td>g7060 </td><td>YM9I_CAEEL^YM9I_CAEEL^Q:165-1547,H:47-510^35.789%ID^E:1.6e-90^RecName: Full=Putative serine protease F56F10.1;^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.</td><td>⋯</td><td>YM9I_CAEEL^YM9I_CAEEL^Q:35-495,H:47-510^35.789%ID^E:7.25e-93^RecName: Full=Putative serine protease F56F10.1;^Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Rhabditina; Rhabditomorpha; Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>PF05577.12^Peptidase_S28^Serine carboxypeptidase S28^52-486^E:8e-146                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         </td><td>sigP:1^21^0.762^YES</td><td>.                                                                                                    </td><td>ENOG410XSGG^protease, serine, 16 (thymus)</td><td>KEGG:cel:CELE_F56F10.1   </td><td>GO:0045121^cellular_component^membrane raft`GO:0008239^molecular_function^dipeptidyl-peptidase activity`GO:0008236^molecular_function^serine-type peptidase activity`GO:0045087^biological_process^innate immune response`GO:0006508^biological_process^proteolysis                                                                                                                                                                                                                                                                                                                                                                                                </td><td>GO:0008236^molecular_function^serine-type peptidase activity`GO:0006508^biological_process^proteolysis                                                                                                                    </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>-3.424802</td><td>1.0781429</td><td>-3.134556</td><td>1.721142e-03</td><td>1.958089e-02</td><td>g7075 </td><td>TSN9_DANRE^TSN9_DANRE^Q:157-810,H:8-230^23.556%ID^E:7.86e-12^RecName: Full=Tetraspanin-9;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii; Teleostei; Ostariophysi; Cypriniformes; Cyprinidae; Danio                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>.</td><td>⋯</td><td>TSN9_DANRE^TSN9_DANRE^Q:12-229,H:8-230^26.222%ID^E:3.91e-25^RecName: Full=Tetraspanin-9;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii; Teleostei; Ostariophysi; Cypriniformes; Cyprinidae; Danio                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>PF00335.20^Tetraspanin^Tetraspanin family^13-226^E:3.8e-43                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   </td><td>.                  </td><td>ExpAA=90.13^PredHel=4^Topology=i17-39o59-81i88-110o203-225i                                          </td><td>ENOG4111IRY^tetraspanin                  </td><td>KEGG:dre:431733`KO:K17350</td><td>GO:0005887^cellular_component^integral component of plasma membrane                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                </td><td>GO:0016021^cellular_component^integral component of membrane                                                                                                                                                              </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>-5.434773</td><td>1.2391311</td><td>-4.243708</td><td>2.198567e-05</td><td>4.265253e-04</td><td>g93314</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.</td><td>⋯</td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            </td><td>.                  </td><td>.                                                                                                    </td><td>.                                        </td><td>.                        </td><td>.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  </td><td>.                                                                                                                                                                                                                         </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10082_length_2005_cov_3.57487_g6803_i1 </td><td> 13.210637</td><td>-2.115922</td><td>0.6939908</td><td>-3.042901</td><td>2.343094e-03</td><td>2.522919e-02</td><td>g6803 </td><td>S35F6_HUMAN^S35F6_HUMAN^Q:12-980,H:18-336^48.457%ID^E:1.23e-96^RecName: Full=Solute carrier family 35 member F6;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </td><td>.</td><td>⋯</td><td>S35F6_HUMAN^S35F6_HUMAN^Q:4-326,H:18-336^48.457%ID^E:5.73e-106^RecName: Full=Solute carrier family 35 member F6;^Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   </td><td>PF08449.11^UAA^UAA transporter family^34-325^E:1e-08`PF08627.10^CRT-like^CRT-like, chloroquine-resistance transporter-like^40-227^E:1.1e-09`PF06027.12^SLC35F^Solute carrier family 35^62-214^E:5.6e-11`PF00892.20^EamA^EamA-like transporter family^85-146^E:5.3e-08                                                                                                                                                                                                                                                                                                                                                                                        </td><td>.                  </td><td>ExpAA=184.31^PredHel=9^Topology=o33-55i75-97o107-124i131-153o163-185i205-224o251-270i291-305o310-332i</td><td>COG0697^membrane                         </td><td>KEGG:hsa:54978           </td><td>GO:0005829^cellular_component^cytosol`GO:0070062^cellular_component^extracellular exosome`GO:0016021^cellular_component^integral component of membrane`GO:0043231^cellular_component^intracellular membrane-bounded organelle`GO:0005765^cellular_component^lysosomal membrane`GO:0005739^cellular_component^mitochondrion`GO:0005654^cellular_component^nucleoplasm`GO:0022857^molecular_function^transmembrane transporter activity`GO:1901029^biological_process^negative regulation of mitochondrial outer membrane permeabilization involved in apoptotic signaling pathway`GO:0008284^biological_process^positive regulation of cell population proliferation</td><td>GO:0055085^biological_process^transmembrane transport`GO:0022857^molecular_function^transmembrane transporter activity`GO:0016021^cellular_component^integral component of membrane`GO:0016020^cellular_component^membrane</td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10087_length_2004_cov_52.6834_g7112_i0 </td><td>  1.431373</td><td>-4.018187</td><td>1.3207673</td><td>-3.009700</td><td>2.615061e-03</td><td>2.759244e-02</td><td>g7112 </td><td>FAS2_SCHAM^FAS2_SCHAM^Q:76-1557,H:8-555^28.905%ID^E:3.41e-63^RecName: Full=Fasciclin-2;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Polyneoptera; Orthoptera; Caelifera; Acrididea; Acridomorpha; Acridoidea; Acrididae; Cyrtacanthacridinae; Schistocerca                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>.</td><td>⋯</td><td>FAS2_SCHAM^FAS2_SCHAM^Q:15-508,H:8-555^28.905%ID^E:6.64e-64^RecName: Full=Fasciclin-2;^Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Polyneoptera; Orthoptera; Caelifera; Acrididea; Acridomorpha; Acridoidea; Acrididae; Cyrtacanthacridinae; Schistocerca                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              </td><td>PF07679.16^I-set^Immunoglobulin I-set domain^138-208^E:4.7e-06`PF13927.6^Ig_3^Immunoglobulin domain^285-368^E:5.6e-14`PF13895.6^Ig_2^Immunoglobulin domain^286-381^E:1.7e-07`PF07679.16^I-set^Immunoglobulin I-set domain^291-381^E:8.7e-14`PF00047.25^ig^Immunoglobulin domain^291-375^E:1e-09`PF13927.6^Ig_3^Immunoglobulin domain^395-466^E:2.7e-09`PF07679.16^I-set^Immunoglobulin I-set domain^395-476^E:1.2e-07`PF00047.25^ig^Immunoglobulin domain^397-475^E:1.5e-07`PF00057.18^Ldl_recept_a^Low-density lipoprotein receptor domain class A^497-532^E:5.9e-12`PF00057.18^Ldl_recept_a^Low-density lipoprotein receptor domain class A^537-572^E:1e-11</td><td>sigP:1^25^0.734^YES</td><td>.                                                                                                    </td><td>.                                        </td><td>.                        </td><td>GO:0016021^cellular_component^integral component of membrane`GO:0007155^biological_process^cell adhesion`GO:0030154^biological_process^cell differentiation`GO:0007399^biological_process^nervous system development                                                                                                                                                                                                                                                                                                                                                                                                                                               </td><td>GO:0005515^molecular_function^protein binding                                                                                                                                                                             </td><td>.</td><td>.</td></tr>
	<tr><td>NODE_10088_length_2004_cov_27.3961_g7113_i0 </td><td> 61.316935</td><td>-7.129926</td><td>1.0023331</td><td>-6.850182</td><td>7.375626e-12</td><td>4.919102e-10</td><td>g7113 </td><td>SANT_PLAFF^SANT_PLAFF^Q:33-527,H:82-243^28.485%ID^E:8.91e-16^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:132-605,H:82-239^30.38%ID^E:9.7e-16^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:66-539,H:82-239^29.747%ID^E:3.17e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:165-638,H:82-239^29.747%ID^E:3.17e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:6-473,H:84-239^28.846%ID^E:3.17e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:99-572,H:82-239^29.747%ID^E:3.54e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:3-407,H:105-239^28.889%ID^E:2.64e-12^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)</td><td>.</td><td>⋯</td><td>SANT_PLAFF^SANT_PLAFF^Q:44-201,H:82-239^30.38%ID^E:4.33e-16^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:11-175,H:82-243^28.485%ID^E:6.49e-16^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:55-212,H:82-239^29.747%ID^E:1.03e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:22-179,H:82-239^29.747%ID^E:1.34e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:33-190,H:82-239^29.747%ID^E:1.5e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:2-157,H:84-239^28.846%ID^E:1.95e-15^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)`SANT_PLAFF^SANT_PLAFF^Q:1-135,H:105-239^28.889%ID^E:1.28e-12^RecName: Full=S-antigen protein;^Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Haemosporida; Plasmodiidae; Plasmodium; Plasmodium (Laverania)</td><td>PF03318.13^ETX_MTX2^Clostridium epsilon toxin ETX/Bacillus mosquitocidal toxin MTX2^366-480^E:0.00011                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        </td><td>.                  </td><td>.                                                                                                    </td><td>.                                        </td><td>.                        </td><td>GO:0020003^cellular_component^symbiont-containing vacuole                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          </td><td>.                                                                                                                                                                                                                         </td><td>.</td><td>.</td></tr>
</tbody>
</table>



 <h2 style="color: #3A9BDC;"> Determine tissue-specific expression  </h2>
 
To identify tissue-specific differential expression (i.e., significantly upregulated genes that are uniquely expressed), each tissue was compared to the other two and tissue-specific genes were extracted from the intersection of the Venn diagram (in Section 9.4.4). For example, the expression in the bioluminescent upper lip was determined by comparing it to both the compound eye and gut tissues.

 <h3 style="color: #3A9BDC;"> Bioluminescent upper lip  </h3>


```R
#create a dataframe with all significantly upregulated genes of the bioluminescent upper lip. 
#merge dataframes that have significantly upregulated genes of the bioluminescent upper lip from pairwise comparisons - Bio Upper Lip vs Gut and Bio Upper Lip vs Eye. 


head(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut) 

head(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye)


```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0</td><td>  7.481691</td><td>4.244836</td><td>1.2023198</td><td>3.570275</td><td>3.566072e-04</td><td>0.0130704556</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1</td><td>144.942119</td><td>3.696833</td><td>0.8140852</td><td>4.529183</td><td>5.921214e-06</td><td>0.0004848792</td></tr>
	<tr><td>NODE_10075_length_2006_cov_50.2255_g7102_i0</td><td> 47.541214</td><td>2.107924</td><td>0.6838391</td><td>3.081562</td><td>2.059176e-03</td><td>0.0439239389</td></tr>
	<tr><td>NODE_10110_length_2001_cov_7.02312_g7128_i0</td><td> 10.296926</td><td>4.431308</td><td>1.0221965</td><td>4.394701</td><td>1.109251e-05</td><td>0.0008282910</td></tr>
	<tr><td>NODE_10314_length_1973_cov_79.1867_g7270_i0</td><td> 51.504733</td><td>1.700903</td><td>0.5199623</td><td>3.269344</td><td>1.077972e-03</td><td>0.0286489771</td></tr>
	<tr><td>NODE_1031_length_4873_cov_1.60834_g724_i0  </td><td> 14.321316</td><td>4.118464</td><td>1.0467162</td><td>3.949132</td><td>7.843501e-05</td><td>0.0040845539</td></tr>
</tbody>
</table>




```R
colnames(sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)[1] <- "gene" 
```


```R
colnames(sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye)[1] <- "gene"
```


```R
#merge
merged_upper_lips_df <- rbind(
  sigOE_UPREGULATED_logfold_Bio_Upper_Lip_Vs_Gut,
  sigOE_UPREGULATED_logfold_Bio_UpperLip_Vs_Eye
)
```


```R
head(merged_upper_lips_df)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td></tr>
</tbody>
</table>




```R
# The same gene can be found in both Bio Upper Lip vs Gut and Bio Upper Lip vs Eye. Remove gene duplicates while retaining one duplicate. 
merged_upper_lips_unique <- merged_upper_lips_df[!duplicated(merged_upper_lips_df$gene), ]


```


```R
head(merged_upper_lips_unique)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td> 32.497969</td><td>4.006190</td><td>1.2628737</td><td>3.190246</td><td>1.421519e-03</td><td>1.671594e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td> 28.846375</td><td>4.198382</td><td>0.7685169</td><td>5.421709</td><td>5.903195e-08</td><td>1.951548e-06</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td>  3.923975</td><td>4.338958</td><td>1.1920667</td><td>3.554069</td><td>3.793190e-04</td><td>5.366714e-03</td></tr>
</tbody>
</table>



<h3 style="color: #3A9BDC;"> Compound eye  </h3>


```R
#compound eye
#create a dataframe with all significantly upregulated genes of the compound eye.  
#merge dataframes that have significantly upregulated genes of the compound eye from pairwise comparisons - Bio Upper Lip vs Eye and Gut vs Eye. 


head(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye)

head(sigOE_DOWNREGULATED_logfold_Gut_vs_Eye)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td> 15.963167</td><td>-3.442293</td><td>0.9226241</td><td> -3.704888</td><td>2.114840e-04</td><td>8.941592e-03</td></tr>
	<tr><td>NODE_101633_length_138_cov_2.06024_g94413_i0</td><td>  7.667808</td><td>-6.149554</td><td>1.0781402</td><td> -5.516364</td><td>3.460855e-08</td><td>4.918433e-06</td></tr>
	<tr><td>NODE_10218_length_1985_cov_18.5554_g7198_i0 </td><td> 61.993605</td><td>-6.085923</td><td>0.6771456</td><td> -8.841580</td><td>9.437437e-19</td><td>3.922413e-16</td></tr>
	<tr><td>NODE_1021_length_4884_cov_10.72_g718_i0     </td><td> 45.575618</td><td>-6.514480</td><td>0.8450182</td><td> -7.517323</td><td>5.590900e-14</td><td>1.734596e-11</td></tr>
	<tr><td>NODE_10232_length_1983_cov_42.0814_g7210_i0 </td><td>434.539937</td><td>-9.286811</td><td>0.7189709</td><td>-12.613262</td><td>1.784501e-36</td><td>6.551496e-33</td></tr>
	<tr><td>NODE_102471_length_136_cov_2.41975_g95251_i0</td><td>  3.353358</td><td>-4.786748</td><td>1.3410759</td><td> -3.571766</td><td>3.545819e-04</td><td>1.301788e-02</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10030_length_2012_cov_88.44_g145_i5    </td><td>77.516326</td><td>-2.004352</td><td>0.7483891</td><td>-2.677266</td><td>7.422565e-03</td><td>4.264394e-02</td></tr>
	<tr><td>NODE_1009_length_4897_cov_2.69228_g709_i0   </td><td>18.842542</td><td>-1.106516</td><td>0.4161180</td><td>-2.659735</td><td>7.820226e-03</td><td>4.417345e-02</td></tr>
	<tr><td>NODE_10106_length_2001_cov_20.3505_g7125_i0 </td><td>32.497969</td><td>-3.786303</td><td>1.2661859</td><td>-3.012530</td><td>2.590799e-03</td><td>1.957288e-02</td></tr>
	<tr><td>NODE_10108_length_2001_cov_9.98304_g7126_i0 </td><td>28.846375</td><td>-4.770626</td><td>0.7722330</td><td>-6.123648</td><td>9.145690e-10</td><td>5.682196e-08</td></tr>
	<tr><td>NODE_101175_length_138_cov_7.48193_g93955_i0</td><td> 3.923975</td><td>-3.824371</td><td>1.2262699</td><td>-3.075752</td><td>2.099721e-03</td><td>1.670540e-02</td></tr>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td>15.963167</td><td>-5.219188</td><td>0.9667613</td><td>-5.307557</td><td>1.111044e-07</td><td>3.921827e-06</td></tr>
</tbody>
</table>




```R
colnames(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye)[1]<- "gene"
```


```R
#merged eye
merged_Eye_df <- rbind(sigOE_DOWNREGULATED_logfold_Bio_UpperLip_Vs_Eye , sigOE_DOWNREGULATED_logfold_Gut_vs_Eye)

```


```R
#The same gene can be found in both Bio Upper Lip vs Eye and Gut vs Eye. Remove gene duplicates while retaining one duplicate. 
merged_Eye_unique <-merged_Eye_df[!duplicated(merged_Eye_df$gene), ]


```


```R
head(merged_Eye_unique)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10126_length_1999_cov_2.60802_g7138_i0 </td><td> 15.963167</td><td>-3.442293</td><td>0.9226241</td><td> -3.704888</td><td>2.114840e-04</td><td>8.941592e-03</td></tr>
	<tr><td>NODE_101633_length_138_cov_2.06024_g94413_i0</td><td>  7.667808</td><td>-6.149554</td><td>1.0781402</td><td> -5.516364</td><td>3.460855e-08</td><td>4.918433e-06</td></tr>
	<tr><td>NODE_10218_length_1985_cov_18.5554_g7198_i0 </td><td> 61.993605</td><td>-6.085923</td><td>0.6771456</td><td> -8.841580</td><td>9.437437e-19</td><td>3.922413e-16</td></tr>
	<tr><td>NODE_1021_length_4884_cov_10.72_g718_i0     </td><td> 45.575618</td><td>-6.514480</td><td>0.8450182</td><td> -7.517323</td><td>5.590900e-14</td><td>1.734596e-11</td></tr>
	<tr><td>NODE_10232_length_1983_cov_42.0814_g7210_i0 </td><td>434.539937</td><td>-9.286811</td><td>0.7189709</td><td>-12.613262</td><td>1.784501e-36</td><td>6.551496e-33</td></tr>
	<tr><td>NODE_102471_length_136_cov_2.41975_g95251_i0</td><td>  3.353358</td><td>-4.786748</td><td>1.3410759</td><td> -3.571766</td><td>3.545819e-04</td><td>1.301788e-02</td></tr>
</tbody>
</table>



<h3 style="color: #3A9BDC;"> Gut  </h3>


```R
#create a dataframe with all significantly upregulated genes of the gut.  
#merge dataframes that have significantly upregulated genes of the gut from pairwise comparisons - Gut vs Eye and Bio Upper Lip vs Gut. 

head(sigOE_UPREGULATED_logfold_Gut_vs_Eye)

head(sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)

```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0 </td><td>  7.481691</td><td>4.329005</td><td>1.1986723</td><td>3.646432</td><td>2.659066e-04</td><td>3.236611e-03</td></tr>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>7.152604</td><td>0.9557768</td><td>7.404069</td><td>1.320736e-13</td><td>2.149895e-11</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>5.960708</td><td>1.2577313</td><td>4.683327</td><td>2.822558e-06</td><td>6.652356e-05</td></tr>
	<tr><td>NODE_10047_length_2010_cov_1.26305_g7081_i0 </td><td>  2.108549</td><td>3.092297</td><td>1.1161236</td><td>2.775470</td><td>5.512192e-03</td><td>3.450172e-02</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>6.415144</td><td>1.3379175</td><td>4.675569</td><td>2.931403e-06</td><td>6.869105e-05</td></tr>
	<tr><td>NODE_1006_length_4905_cov_1.0567_g707_i0    </td><td>  2.237450</td><td>3.482954</td><td>1.2248812</td><td>2.814833</td><td>4.880265e-03</td><td>3.151585e-02</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>-6.417432</td><td>0.9300257</td><td>-6.864408</td><td>6.676754e-12</td><td>4.486477e-10</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>-3.424802</td><td>1.0781429</td><td>-3.134556</td><td>1.721142e-03</td><td>1.958089e-02</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>-5.434773</td><td>1.2391311</td><td>-4.243708</td><td>2.198567e-05</td><td>4.265253e-04</td></tr>
	<tr><td>NODE_10082_length_2005_cov_3.57487_g6803_i1 </td><td> 13.210637</td><td>-2.115922</td><td>0.6939908</td><td>-3.042901</td><td>2.343094e-03</td><td>2.522919e-02</td></tr>
	<tr><td>NODE_10087_length_2004_cov_52.6834_g7112_i0 </td><td>  1.431373</td><td>-4.018187</td><td>1.3207673</td><td>-3.009700</td><td>2.615061e-03</td><td>2.759244e-02</td></tr>
	<tr><td>NODE_10088_length_2004_cov_27.3961_g7113_i0 </td><td> 61.316935</td><td>-7.129926</td><td>1.0023331</td><td>-6.850182</td><td>7.375626e-12</td><td>4.919102e-10</td></tr>
</tbody>
</table>




```R
colnames(sigOE_UPREGULATED_logfold_Gut_vs_Eye)[1] <- "gene"
```


```R
colnames(sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)[1]<- "gene"
```


```R
#merge
merged_Gut_df <- rbind(sigOE_UPREGULATED_logfold_Gut_vs_Eye , sigOE_DOWNREGULATED_logfold_Bio_Upper_Lip_Vs_Gut)

#The same gene can be found in both Gut vs Eye and Bio Upper Lip vs Eye. Remove gene duplicates while retaining one duplicate.
merged_Gut_unique <-merged_Gut_df[!duplicated(merged_Gut_df$gene), ]


```


```R
head(merged_Gut_unique)
```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>gene</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10010_length_2015_cov_7.32245_g7057_i0 </td><td>  7.481691</td><td>4.329005</td><td>1.1986723</td><td>3.646432</td><td>2.659066e-04</td><td>3.236611e-03</td></tr>
	<tr><td>NODE_10015_length_2014_cov_682.243_g7060_i0 </td><td>570.345545</td><td>7.152604</td><td>0.9557768</td><td>7.404069</td><td>1.320736e-13</td><td>2.149895e-11</td></tr>
	<tr><td>NODE_10040_length_2011_cov_2.41973_g7075_i0 </td><td> 14.886486</td><td>5.960708</td><td>1.2577313</td><td>4.683327</td><td>2.822558e-06</td><td>6.652356e-05</td></tr>
	<tr><td>NODE_10047_length_2010_cov_1.26305_g7081_i0 </td><td>  2.108549</td><td>3.092297</td><td>1.1161236</td><td>2.775470</td><td>5.512192e-03</td><td>3.450172e-02</td></tr>
	<tr><td>NODE_100534_length_140_cov_1.47059_g93314_i0</td><td> 21.181026</td><td>6.415144</td><td>1.3379175</td><td>4.675569</td><td>2.931403e-06</td><td>6.869105e-05</td></tr>
	<tr><td>NODE_1006_length_4905_cov_1.0567_g707_i0    </td><td>  2.237450</td><td>3.482954</td><td>1.2248812</td><td>2.814833</td><td>4.880265e-03</td><td>3.151585e-02</td></tr>
</tbody>
</table>



<h3 style="color: #3A9BDC;"> Venn Diagram - Extract tissue-specific genes   </h3>


```R
#generate a venn diagram to visualize shared significantly upregulated genes across tissue types and extract the genes that are unique to each tissue type. 
unique_venn_list <- list(
  Bio_Upper_Lip = merged_upper_lips_unique$gene  , 
  Gut = merged_Gut_unique$gene,
  Compound_Eye = merged_Eye_unique$gene
)

ggvenn_unique <- ggvenn(
  unique_venn_list, 
  fill_color = c('#AC97C9','#C97D97', '#F2C93D'),
  stroke_size = .7, set_name_size = 6, text_size = 5
)

ggvenn_unique 
```


![png](output_231_0.png)



```R
# Open a PDF device
pdf("Luminous_DGE.pdf", width = 8, height = 6)


ggvenn_unique 

dev.off()

```


<strong>png:</strong> 2



```R
#prep dataframes for extraction
Bio_Upper_Lip <- as.data.frame(merged_upper_lips_unique$gene)
colnames(Bio_Upper_Lip)[1]<- "gene"
Gut <- as.data.frame(merged_Gut_unique$gene)
colnames(Gut)[1]<- "gene"
Compound_eye <- as.data.frame(merged_Eye_unique$gene)
colnames(Compound_eye)[1]<- "gene"
```


```R
# compare and extract unique genes for each tissue type
unique_genes_bio_upper_lip <- anti_join(Bio_Upper_Lip, Gut, by = "gene") %>%
  anti_join(Compound_eye, by = "gene")

unique_genes_gut <- anti_join(Gut, Bio_Upper_Lip, by = "gene") %>%
  anti_join(Compound_eye, by = "gene")

unique_genes_compound_eye <- anti_join(Compound_eye, Bio_Upper_Lip, by = "gene") %>%
  anti_join(Gut, by = "gene")


```


```R
nrow(unique_genes_bio_upper_lip)
```


595



```R
nrow(unique_genes_gut)
```


2534



```R
nrow(unique_genes_compound_eye)
```


1144



```R
#add the annotations back to the unique genes in each tissue type by subsetting 
```


```R
unique_genes_bio_upper_lip_info  <- merged_upper_lips_unique %>%
  filter(gene %in% unique_genes_bio_upper_lip$gene)
```


```R
nrow(unique_genes_bio_upper_lip_info)
```


595



```R
unique_genes_eye_info  <- merged_Eye_unique %>%
  filter(gene %in% unique_genes_compound_eye$gene)
```


```R
nrow(unique_genes_eye_info)
```


1144



```R
unique_genes_gut_info  <- merged_Gut_unique %>%
  filter(gene %in% unique_genes_gut$gene)
```


```R
nrow(unique_genes_gut_info)
```


2534



```R
#add annotations back 
colnames(unique_genes_bio_upper_lip_info)[1]<- "transcript_id"
unique_genes_bio_upper_lip_info_annot <- left_join(unique_genes_bio_upper_lip_info,Trinotate_lym_subset,by="transcript_id") 

```


```R
#add annotations back 
colnames(unique_genes_eye_info)[1]<- "transcript_id"
unique_genes_eye_info_annot <- left_join(unique_genes_eye_info,Trinotate_lym_subset,by="transcript_id") 
```


```R
#add annotations backs
colnames(unique_genes_gut_info)[1] <- "transcript_id"
unique_genes_gut_info_annot <- left_join(unique_genes_gut_info,Trinotate_lym_subset,by="transcript_id") 

```


```R
#write.csv(unique_genes_bio_upper_lip_info_annot, file = "df_Vtsujii_sigfig_upreg_unique_BioUpperLip.csv")
```


```R
#write.csv(unique_genes_eye_info_annot, file = "df_Vtsujii_sigfig_upreg_unique_comEye.csv")
```


```R
#write.csv(unique_genes_gut_info_annot, file = "df_Vtsujii_sigfig_upreg_unique_Gut.csv")
```

 <h1 style="color: #3A9BDC;"> DGE - GO enrichment analyses for bioluminescent upper lip  </h1> 
 
The topGO package was used to perform GO enrichment analysis on significantly upregulated genes (i.e., uniquely expressed) in the bioluminescent upper lip. (Alexa Rahnenfuhrer, 2024). The figures representing the GO enrichments presented here were utilized for analysis, but were not included in the main text of the publication. The GO enrichments were visualized with GO-Figure! package for publication (Reijnders and Waterhouse, 2021). GO-Figure! reference https://github.com/lmesrop/BCN_publication/tree/main/Go-Figure!. 



```R
#import the trinotate go sheet from Trinotate output 
geneID2GO_bio <- readMappings(file ="Trinotate_go_lym.txt")

```


```R
#significantly upregulated genes (i.e. expressed uniquely) in the bioluminescent upper lip 

head(unique_genes_bio_upper_lip_info)

```


<table class="dataframe">
<caption>A tibble: 6 × 7</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>NODE_10035_length_2011_cov_7.68814_g7072_i0 </td><td> 31.497464</td><td>3.549189</td><td>1.1116132</td><td>3.185621</td><td>1.444434e-03</td><td>1.694083e-02</td></tr>
	<tr><td>NODE_10049_length_2009_cov_1010.67_g4245_i1 </td><td>144.942119</td><td>3.780206</td><td>0.8027568</td><td>4.702933</td><td>2.564508e-06</td><td>6.123575e-05</td></tr>
	<tr><td>NODE_10092_length_2003_cov_126.275_g7115_i0 </td><td> 87.910226</td><td>4.226131</td><td>0.9071064</td><td>4.653001</td><td>3.271381e-06</td><td>7.653491e-05</td></tr>
	<tr><td>NODE_10123_length_1999_cov_24.5818_g6951_i1 </td><td> 12.970945</td><td>1.921576</td><td>0.6823847</td><td>2.815052</td><td>4.876929e-03</td><td>4.595970e-02</td></tr>
	<tr><td>NODE_10147_length_1995_cov_11.3794_g7152_i0 </td><td>  7.942492</td><td>3.961235</td><td>1.1603461</td><td>3.372871</td><td>7.438888e-04</td><td>9.602505e-03</td></tr>
	<tr><td>NODE_102011_length_137_cov_2.19512_g94791_i0</td><td> 83.590549</td><td>6.076569</td><td>1.6982637</td><td>3.942479</td><td>8.064363e-05</td><td>1.377157e-03</td></tr>
</tbody>
</table>




```R
#save the transcript ids of all the annotated genes under geneNames object 
geneNames_bio<- as.character(Trinotate_lym_subset$transcript_id)
```


```R
#save the transcript 
myInterestingGenes_bio= as.character(unique_genes_bio_upper_lip_info$gene)

```


```R
#subset the genesNames by the transcript IDs in my red module 
geneList_bio <- factor(as.integer(geneNames_bio %in% myInterestingGenes_bio))
names(geneList_bio) <- geneNames_bio
head(geneList_bio)

```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>NODE_100000_length_141_cov_1.05814_g92780_i0</dt><dd>0</dd><dt>NODE_100001_length_141_cov_1.03488_g92781_i0</dt><dd>0</dd><dt>NODE_100002_length_141_cov_0.325581_g92782_i0</dt><dd>0</dd><dt>NODE_100003_length_141_cov_0.0116279_g92783_i0</dt><dd>0</dd><dt>NODE_100004_length_141_cov_0_g92784_i0</dt><dd>0</dd><dt>NODE_100005_length_141_cov_0_g92785_i0</dt><dd>0</dd></dl>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'0'</li><li>'1'</li></ol>
</details>


<h2 style="color: #3A9BDC;"> Use topGO to identify enriched biological processes in the bioluminescent upper lip  </h2> 


```R
#run the topGO function.
GOdata_bio <- new("topGOdata", ontology = "BP", allGenes = geneList_bio,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO_bio)

```

    
    Building most specific GOs .....
    
    	( 12645 GO terms found. )
    
    
    Build GO DAG topology ..........
    
    	( 13427 GO terms and 31130 relations. )
    
    
    Annotating nodes ...............
    
    	( 15506 genes annotated to the GO terms. )
    



```R
results_go_bio <- runTest(GOdata_bio, algorithm="weight01", statistic="Fisher")
```

    
    			 -- Weight01 Algorithm -- 
    
    		 the algorithm is scoring 1881 nontrivial nodes
    		 parameters: 
    			 test statistic: fisher
    
    
    	 Level 16:	1 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 15:	4 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 14:	8 nodes to be scored	(6 eliminated genes)
    
    
    	 Level 13:	25 nodes to be scored	(115 eliminated genes)
    
    
    	 Level 12:	56 nodes to be scored	(265 eliminated genes)
    
    
    	 Level 11:	94 nodes to be scored	(1648 eliminated genes)
    
    
    	 Level 10:	139 nodes to be scored	(3396 eliminated genes)
    
    
    	 Level 9:	203 nodes to be scored	(4708 eliminated genes)
    
    
    	 Level 8:	243 nodes to be scored	(5830 eliminated genes)
    
    
    	 Level 7:	302 nodes to be scored	(8018 eliminated genes)
    
    
    	 Level 6:	315 nodes to be scored	(11122 eliminated genes)
    
    
    	 Level 5:	254 nodes to be scored	(12818 eliminated genes)
    
    
    	 Level 4:	149 nodes to be scored	(14217 eliminated genes)
    
    
    	 Level 3:	70 nodes to be scored	(14877 eliminated genes)
    
    
    	 Level 2:	17 nodes to be scored	(15171 eliminated genes)
    
    
    	 Level 1:	1 nodes to be scored	(15489 eliminated genes)
    



```R
#retrieve the GO enrichment 
goEnrichment_bio   <- GenTable(GOdata_bio, Fisher = results_go_bio, orderBy = "Fisher", topNodes =100, numChar =1000 )
```


```R
goEnrichment_bio$Fisher <- as.numeric(goEnrichment_bio$Fisher)
goEnrichment_bio <- goEnrichment_bio[goEnrichment_bio$Fisher < 0.05,] 
goEnrichment_bio <- goEnrichment_bio[goEnrichment_bio$Significant > 1,] 
#goEnrichment_bio <- goEnrichment_bio[,c("GO.ID","Term", "Annotated","Significant", "Expected", "Fisher")]

```


```R
goEnrichment_bio
```


<table class="dataframe">
<caption>A data.frame: 48 × 6</caption>
<thead>
	<tr><th></th><th scope=col>GO.ID</th><th scope=col>Term</th><th scope=col>Annotated</th><th scope=col>Significant</th><th scope=col>Expected</th><th scope=col>Fisher</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0006508</td><td>proteolysis                                             </td><td>1540</td><td>25</td><td>13.41</td><td>3.400e-09</td></tr>
	<tr><th scope=row>2</th><td>GO:0008218</td><td>bioluminescence                                         </td><td>  46</td><td> 6</td><td> 0.40</td><td>2.700e-06</td></tr>
	<tr><th scope=row>3</th><td>GO:0032504</td><td>multicellular organism reproduction                     </td><td> 969</td><td>10</td><td> 8.44</td><td>8.700e-06</td></tr>
	<tr><th scope=row>4</th><td>GO:0001519</td><td>peptide amidation                                       </td><td>  19</td><td> 4</td><td> 0.17</td><td>1.900e-05</td></tr>
	<tr><th scope=row>5</th><td>GO:0006030</td><td>chitin metabolic process                                </td><td>  76</td><td> 6</td><td> 0.66</td><td>5.200e-05</td></tr>
	<tr><th scope=row>6</th><td>GO:0018401</td><td>peptidyl-proline hydroxylation to 4-hydroxy-L-proline   </td><td>  25</td><td> 4</td><td> 0.22</td><td>6.000e-05</td></tr>
	<tr><th scope=row>7</th><td>GO:0032849</td><td>positive regulation of cellular pH reduction            </td><td>  10</td><td> 3</td><td> 0.09</td><td>7.400e-05</td></tr>
	<tr><th scope=row>8</th><td>GO:0006556</td><td>S-adenosylmethionine biosynthetic process               </td><td>  10</td><td> 3</td><td> 0.09</td><td>7.400e-05</td></tr>
	<tr><th scope=row>9</th><td>GO:2001225</td><td>regulation of chloride transport                        </td><td>  10</td><td> 3</td><td> 0.09</td><td>7.400e-05</td></tr>
	<tr><th scope=row>10</th><td>GO:0009620</td><td>response to fungus                                      </td><td>  60</td><td> 5</td><td> 0.52</td><td>9.400e-05</td></tr>
	<tr><th scope=row>11</th><td>GO:0032230</td><td>positive regulation of synaptic transmission, GABAergic </td><td>  11</td><td> 3</td><td> 0.10</td><td>1.000e-04</td></tr>
	<tr><th scope=row>12</th><td>GO:0007586</td><td>digestion                                               </td><td>  93</td><td> 7</td><td> 0.81</td><td>1.100e-04</td></tr>
	<tr><th scope=row>13</th><td>GO:0036378</td><td>calcitriol biosynthetic process from calciol            </td><td>   3</td><td> 2</td><td> 0.03</td><td>2.200e-04</td></tr>
	<tr><th scope=row>14</th><td>GO:0010164</td><td>response to cesium ion                                  </td><td>   3</td><td> 2</td><td> 0.03</td><td>2.200e-04</td></tr>
	<tr><th scope=row>15</th><td>GO:0048771</td><td>tissue remodeling                                       </td><td>  76</td><td> 6</td><td> 0.66</td><td>3.100e-04</td></tr>
	<tr><th scope=row>16</th><td>GO:0042730</td><td>fibrinolysis                                            </td><td>  17</td><td> 3</td><td> 0.15</td><td>4.000e-04</td></tr>
	<tr><th scope=row>17</th><td>GO:2001150</td><td>positive regulation of dipeptide transmembrane transport</td><td>   4</td><td> 2</td><td> 0.03</td><td>4.500e-04</td></tr>
	<tr><th scope=row>18</th><td>GO:0015670</td><td>carbon dioxide transport                                </td><td>   5</td><td> 2</td><td> 0.04</td><td>7.400e-04</td></tr>
	<tr><th scope=row>19</th><td>GO:0044719</td><td>regulation of imaginal disc-derived wing size           </td><td>  49</td><td> 4</td><td> 0.43</td><td>8.600e-04</td></tr>
	<tr><th scope=row>20</th><td>GO:0038166</td><td>angiotensin-activated signaling pathway                 </td><td>   6</td><td> 2</td><td> 0.05</td><td>1.100e-03</td></tr>
	<tr><th scope=row>21</th><td>GO:0045780</td><td>positive regulation of bone resorption                  </td><td>   8</td><td> 2</td><td> 0.07</td><td>2.040e-03</td></tr>
	<tr><th scope=row>22</th><td>GO:0007218</td><td>neuropeptide signaling pathway                          </td><td>  32</td><td> 3</td><td> 0.28</td><td>2.660e-03</td></tr>
	<tr><th scope=row>23</th><td>GO:0010043</td><td>response to zinc ion                                    </td><td>  41</td><td> 3</td><td> 0.36</td><td>5.400e-03</td></tr>
	<tr><th scope=row>24</th><td>GO:0001658</td><td>branching involved in ureteric bud morphogenesis        </td><td>  13</td><td> 2</td><td> 0.11</td><td>5.510e-03</td></tr>
	<tr><th scope=row>25</th><td>GO:0045672</td><td>positive regulation of osteoclast differentiation       </td><td>  14</td><td> 2</td><td> 0.12</td><td>6.390e-03</td></tr>
	<tr><th scope=row>26</th><td>GO:0006730</td><td>one-carbon metabolic process                            </td><td>  45</td><td> 3</td><td> 0.39</td><td>7.010e-03</td></tr>
	<tr><th scope=row>29</th><td>GO:0009405</td><td>pathogenesis                                            </td><td>  20</td><td> 2</td><td> 0.17</td><td>1.290e-02</td></tr>
	<tr><th scope=row>30</th><td>GO:0071498</td><td>cellular response to fluid shear stress                 </td><td>  20</td><td> 2</td><td> 0.17</td><td>1.290e-02</td></tr>
	<tr><th scope=row>31</th><td>GO:0007613</td><td>memory                                                  </td><td> 109</td><td> 4</td><td> 0.95</td><td>1.514e-02</td></tr>
	<tr><th scope=row>32</th><td>GO:0030574</td><td>collagen catabolic process                              </td><td>  22</td><td> 2</td><td> 0.19</td><td>1.551e-02</td></tr>
	<tr><th scope=row>33</th><td>GO:0042475</td><td>odontogenesis of dentin-containing tooth                </td><td>  23</td><td> 2</td><td> 0.20</td><td>1.689e-02</td></tr>
	<tr><th scope=row>37</th><td>GO:0015701</td><td>bicarbonate transport                                   </td><td>  25</td><td> 2</td><td> 0.22</td><td>1.980e-02</td></tr>
	<tr><th scope=row>38</th><td>GO:0007585</td><td>respiratory gaseous exchange                            </td><td>  25</td><td> 2</td><td> 0.22</td><td>1.980e-02</td></tr>
	<tr><th scope=row>39</th><td>GO:0055085</td><td>transmembrane transport                                 </td><td>1050</td><td>15</td><td> 9.14</td><td>2.183e-02</td></tr>
	<tr><th scope=row>40</th><td>GO:0055114</td><td>oxidation-reduction process                             </td><td>1003</td><td>15</td><td> 8.73</td><td>2.223e-02</td></tr>
	<tr><th scope=row>41</th><td>GO:0042738</td><td>exogenous drug catabolic process                        </td><td>  27</td><td> 2</td><td> 0.24</td><td>2.291e-02</td></tr>
	<tr><th scope=row>42</th><td>GO:0015771</td><td>trehalose transport                                     </td><td>  27</td><td> 2</td><td> 0.24</td><td>2.291e-02</td></tr>
	<tr><th scope=row>43</th><td>GO:0046903</td><td>secretion                                               </td><td> 755</td><td> 8</td><td> 6.57</td><td>2.302e-02</td></tr>
	<tr><th scope=row>44</th><td>GO:0003073</td><td>regulation of systemic arterial blood pressure          </td><td>  28</td><td> 2</td><td> 0.24</td><td>2.453e-02</td></tr>
	<tr><th scope=row>45</th><td>GO:0042359</td><td>vitamin D metabolic process                             </td><td>   6</td><td> 3</td><td> 0.05</td><td>2.552e-02</td></tr>
	<tr><th scope=row>54</th><td>GO:0030206</td><td>chondroitin sulfate biosynthetic process                </td><td>  31</td><td> 2</td><td> 0.27</td><td>2.967e-02</td></tr>
	<tr><th scope=row>55</th><td>GO:0010092</td><td>specification of animal organ identity                  </td><td>  21</td><td> 2</td><td> 0.18</td><td>3.416e-02</td></tr>
	<tr><th scope=row>66</th><td>GO:0043434</td><td>response to peptide hormone                             </td><td> 232</td><td> 5</td><td> 2.02</td><td>3.459e-02</td></tr>
	<tr><th scope=row>67</th><td>GO:0009268</td><td>response to pH                                          </td><td>  37</td><td> 2</td><td> 0.32</td><td>4.109e-02</td></tr>
	<tr><th scope=row>68</th><td>GO:0006629</td><td>lipid metabolic process                                 </td><td>1154</td><td> 9</td><td>10.05</td><td>4.271e-02</td></tr>
	<tr><th scope=row>79</th><td>GO:1902017</td><td>regulation of cilium assembly                           </td><td>  38</td><td> 2</td><td> 0.33</td><td>4.313e-02</td></tr>
	<tr><th scope=row>80</th><td>GO:0008543</td><td>fibroblast growth factor receptor signaling pathway     </td><td>  38</td><td> 2</td><td> 0.33</td><td>4.313e-02</td></tr>
	<tr><th scope=row>81</th><td>GO:0002009</td><td>morphogenesis of an epithelium                          </td><td> 558</td><td>13</td><td> 4.86</td><td>4.514e-02</td></tr>
</tbody>
</table>




```R
#write.table(goEnrichment_bio, "df_TopGO_Vargula_tsujii_DE_unique_BioUpperLip_BP.tsv",sep = "\t", quote=FALSE)

```


```R
myterms_bio =goEnrichment_bio$GO.ID 
mygenes_bio = genesInTerm(GOdata_bio, myterms_bio)
```


```R
#extract the transcript ids for each GO term
for (i in 1:length(myterms_bio))
{
   myterm_bio <- myterms_bio[i]
   mygenesforterm_bio <- mygenes_bio[myterm_bio][[1]]
   myfactor_bio <- mygenesforterm_bio %in% myInterestingGenes_bio 
   mygenesforterm2_bio <- mygenesforterm_bio[myfactor_bio == TRUE]
   mygenesforterm2_bio <- paste(mygenesforterm2_bio, collapse=',')
   print(paste("Term",myterm_bio,"genes:",mygenesforterm2_bio))
}

```


```R
ntop = 48
ggdata <- goEnrichment_bio[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) 
plot_Bio_UL_BP <- ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = Significant, fill = -log10(Fisher))) +

 ylim(1,11) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  labs(
    title = 'GO Analysis - Biological Process')+
 theme_bw(base_size = 24) +
labs(size= "Number of Genes")+
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

   axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 13, face = 'bold', vjust = 0.5),
    #axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 8, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), 
    legend.text = element_text(size = 14, face = "bold"), 
    title = element_text(size = 14, face = "bold")) +



  coord_flip()
#dev.off()
```


```R
plot_Bio_UL_BP + labs(x = NULL)
```


![png](output_267_0.png)



```R
options(repr.plot.width=10, repr.plot.height=8, repr.plot.res = 500)
```

 <h2 style="color: #3A9BDC;"> Use topGO to identify enriched molecular functions in the bioluminescent upper lip </h2>


```R
#run the topGO function.
GOdata_bio_MF <- new("topGOdata", ontology = "MF", allGenes = geneList_bio,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO_bio)

```

    
    Building most specific GOs .....
    
    	( 3583 GO terms found. )
    
    
    Build GO DAG topology ..........
    
    	( 3618 GO terms and 4733 relations. )
    
    
    Annotating nodes ...............
    
    	( 16040 genes annotated to the GO terms. )
    



```R
results_go_bio_MF <- runTest(GOdata_bio_MF, algorithm="weight01", statistic="Fisher")
```

    
    			 -- Weight01 Algorithm -- 
    
    		 the algorithm is scoring 401 nontrivial nodes
    		 parameters: 
    			 test statistic: fisher
    
    
    	 Level 12:	1 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 11:	2 nodes to be scored	(0 eliminated genes)
    
    
    	 Level 10:	3 nodes to be scored	(7 eliminated genes)
    
    
    	 Level 9:	11 nodes to be scored	(38 eliminated genes)
    
    
    	 Level 8:	28 nodes to be scored	(58 eliminated genes)
    
    
    	 Level 7:	54 nodes to be scored	(2776 eliminated genes)
    
    
    	 Level 6:	79 nodes to be scored	(3524 eliminated genes)
    
    
    	 Level 5:	92 nodes to be scored	(6256 eliminated genes)
    
    
    	 Level 4:	88 nodes to be scored	(8602 eliminated genes)
    
    
    	 Level 3:	33 nodes to be scored	(13147 eliminated genes)
    
    
    	 Level 2:	9 nodes to be scored	(14316 eliminated genes)
    
    
    	 Level 1:	1 nodes to be scored	(15865 eliminated genes)
    



```R
#retrieve the GO enrichment 
goEnrichment_bio_MF   <- GenTable(GOdata_bio_MF, Fisher = results_go_bio_MF, orderBy = "Fisher", topNodes =100, numChar =1000 )
```


```R
goEnrichment_bio_MF$Fisher <- as.numeric(goEnrichment_bio_MF$Fisher)
goEnrichment_bio_MF <- goEnrichment_bio_MF[goEnrichment_bio_MF$Fisher < 0.05,] 
goEnrichment_bio_MF <- goEnrichment_bio_MF[goEnrichment_bio_MF$Significant > 1,] 
goEnrichment_bio_MF <- goEnrichment_bio_MF[,c("GO.ID","Term", "Annotated","Significant", "Expected", "Fisher")]

```


```R
goEnrichment_bio_MF 


```


<table class="dataframe">
<caption>A data.frame: 29 × 6</caption>
<thead>
	<tr><th></th><th scope=col>GO.ID</th><th scope=col>Term</th><th scope=col>Annotated</th><th scope=col>Significant</th><th scope=col>Expected</th><th scope=col>Fisher</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0004252</td><td>serine-type endopeptidase activity                                                                                           </td><td>241</td><td>13</td><td>2.28</td><td>0.0000005</td></tr>
	<tr><th scope=row>2</th><td>GO:0004222</td><td>metalloendopeptidase activity                                                                                                </td><td>118</td><td> 9</td><td>1.12</td><td>0.0000018</td></tr>
	<tr><th scope=row>3</th><td>GO:0004089</td><td>carbonate dehydratase activity                                                                                               </td><td> 26</td><td> 5</td><td>0.25</td><td>0.0000040</td></tr>
	<tr><th scope=row>4</th><td>GO:0047712</td><td>Cypridina-luciferin 2-monooxygenase activity                                                                                 </td><td> 46</td><td> 6</td><td>0.44</td><td>0.0000045</td></tr>
	<tr><th scope=row>5</th><td>GO:0004478</td><td>methionine adenosyltransferase activity                                                                                      </td><td>  7</td><td> 3</td><td>0.07</td><td>0.0000280</td></tr>
	<tr><th scope=row>6</th><td>GO:0004598</td><td>peptidylamidoglycolate lyase activity                                                                                        </td><td> 20</td><td> 4</td><td>0.19</td><td>0.0000330</td></tr>
	<tr><th scope=row>7</th><td>GO:0004504</td><td>peptidylglycine monooxygenase activity                                                                                       </td><td> 20</td><td> 4</td><td>0.19</td><td>0.0000330</td></tr>
	<tr><th scope=row>8</th><td>GO:0008061</td><td>chitin binding                                                                                                               </td><td>119</td><td> 7</td><td>1.13</td><td>0.0001400</td></tr>
	<tr><th scope=row>9</th><td>GO:0030343</td><td>vitamin D3 25-hydroxylase activity                                                                                           </td><td>  3</td><td> 2</td><td>0.03</td><td>0.0002700</td></tr>
	<tr><th scope=row>10</th><td>GO:0004867</td><td>serine-type endopeptidase inhibitor activity                                                                                 </td><td> 95</td><td> 6</td><td>0.90</td><td>0.0002800</td></tr>
	<tr><th scope=row>11</th><td>GO:0004181</td><td>metallocarboxypeptidase activity                                                                                             </td><td> 36</td><td> 4</td><td>0.34</td><td>0.0003600</td></tr>
	<tr><th scope=row>12</th><td>GO:0016831</td><td>carboxy-lyase activity                                                                                                       </td><td>101</td><td> 6</td><td>0.96</td><td>0.0004000</td></tr>
	<tr><th scope=row>13</th><td>GO:0005184</td><td>neuropeptide hormone activity                                                                                                </td><td>  5</td><td> 3</td><td>0.05</td><td>0.0005200</td></tr>
	<tr><th scope=row>14</th><td>GO:0004574</td><td>oligo-1,6-glucosidase activity                                                                                               </td><td>  5</td><td> 2</td><td>0.05</td><td>0.0008800</td></tr>
	<tr><th scope=row>15</th><td>GO:0004575</td><td>sucrose alpha-glucosidase activity                                                                                           </td><td>  5</td><td> 2</td><td>0.05</td><td>0.0008800</td></tr>
	<tr><th scope=row>16</th><td>GO:0004656</td><td>procollagen-proline 4-dioxygenase activity                                                                                   </td><td> 28</td><td> 3</td><td>0.27</td><td>0.0023000</td></tr>
	<tr><th scope=row>17</th><td>GO:0004064</td><td>arylesterase activity                                                                                                        </td><td>  9</td><td> 2</td><td>0.09</td><td>0.0030700</td></tr>
	<tr><th scope=row>18</th><td>GO:0005507</td><td>copper ion binding                                                                                                           </td><td> 65</td><td> 4</td><td>0.62</td><td>0.0033500</td></tr>
	<tr><th scope=row>19</th><td>GO:0046422</td><td>violaxanthin de-epoxidase activity                                                                                           </td><td> 10</td><td> 2</td><td>0.09</td><td>0.0038200</td></tr>
	<tr><th scope=row>20</th><td>GO:0031418</td><td>L-ascorbic acid binding                                                                                                      </td><td> 41</td><td> 3</td><td>0.39</td><td>0.0068300</td></tr>
	<tr><th scope=row>21</th><td>GO:0005509</td><td>calcium ion binding                                                                                                          </td><td>640</td><td>13</td><td>6.06</td><td>0.0078900</td></tr>
	<tr><th scope=row>23</th><td>GO:0080030</td><td>methyl indole-3-acetate esterase activity                                                                                    </td><td> 17</td><td> 2</td><td>0.16</td><td>0.0110500</td></tr>
	<tr><th scope=row>25</th><td>GO:0015459</td><td>potassium channel regulator activity                                                                                         </td><td> 24</td><td> 2</td><td>0.23</td><td>0.0214800</td></tr>
	<tr><th scope=row>26</th><td>GO:0016702</td><td>oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen</td><td> 63</td><td> 3</td><td>0.60</td><td>0.0219000</td></tr>
	<tr><th scope=row>27</th><td>GO:0015151</td><td>alpha-glucoside transmembrane transporter activity                                                                           </td><td> 26</td><td> 2</td><td>0.25</td><td>0.0249900</td></tr>
	<tr><th scope=row>28</th><td>GO:0015574</td><td>trehalose transmembrane transporter activity                                                                                 </td><td> 26</td><td> 2</td><td>0.25</td><td>0.0249900</td></tr>
	<tr><th scope=row>29</th><td>GO:0018833</td><td>DDT-dehydrochlorinase activity                                                                                               </td><td> 26</td><td> 2</td><td>0.25</td><td>0.0249900</td></tr>
	<tr><th scope=row>34</th><td>GO:0005506</td><td>iron ion binding                                                                                                             </td><td>265</td><td> 6</td><td>2.51</td><td>0.0405900</td></tr>
	<tr><th scope=row>35</th><td>GO:0008395</td><td>steroid hydroxylase activity                                                                                                 </td><td> 34</td><td> 2</td><td>0.32</td><td>0.0410700</td></tr>
</tbody>
</table>




```R
ntop = 29
ggdata <- goEnrichment_bio_MF[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) 
plot_Bio_UL_MF <- ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = Significant, fill = -log10(Fisher))) +

 ylim(1,11) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  labs(
    title = 'GO Analysis - Molecular Function')+
 theme_bw(base_size = 24) +
labs(size= "Number of Genes")+
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

   axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 13, face = 'bold', vjust = 0.5),
    #axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_blank(),
    #axis.title.y = element_text(size = 8, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), 
    legend.text = element_text(size = 14, face = "bold"), 
    title = element_text(size = 14, face = "bold")) +



  coord_flip()
#dev.off()
```


```R
plot_Bio_UL_MF + labs(x=NULL)
```


![png](output_276_0.png)



```R
options(repr.plot.width=16, repr.plot.height=8, repr.plot.res = 500)
```

 <h1 style="color: #3A9BDC;"> References  </h1> 
 
P. Langfelder, S. Horvath, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008)

M. I. Love, W. Huber, S. Anders, Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550 (2014)
 
Alexa A, Rahnenfuhrer J, topGO: Enrichment Analysis for Gene Ontology (2024)

Shannon, Paul et al., Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research vol. 13,11 (2003)


```R

```
