
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: Fit WGCNA 
#
# Description: Fit WGCNA networks to multi-strain A.fumigatus RNA-seq dataset, containing normia & hypoxia samples 
#
# Lab: Cramer 
# Author: Owen WIlkins 
# Edited by: Charles Puerner
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(WGCNA)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(igraph)
ppi=300

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#setwd in working directory.

# read in normalized counts 
norm_counts <- readRDS("./outputData/Exploratory/normalized_counts_vst.rds")

# read in phenotypic data 
meta <- read.csv("inputData/metadata.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
# meta$sample <- NULL
meta$growth <- factor(meta$growth, levels = c("normoxia_planktonic","hypoxia_planktonic","biofilm"))
meta$condition <- factor(meta$condition, levels = c("normoxia","hypoxia"))
meta$batch_no <- factor(as.character(meta$batch_no), levels = c("7",  "2",  "3",  "4",  "1",  "B1", "B2", "B3", "B4"))
meta$group <- factor(meta$group, levels = unique(meta$group))
meta$culture <- factor(meta$culture, levels = (c("planktonic","biofilm")))

# remove samples failing QC 
meta <- meta[!meta$qc == "FAIL",]

# subset counts for samples passing QC 
norm_counts <- norm_counts[rownames(norm_counts) %in% rownames(meta),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# determine soft-thresholding power for WGCNA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run pick soft threshold 
sft <- pickSoftThreshold(norm_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed")

# set constants 
cex1 = 0.9
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

dir.create("./outputData/WGCNA/plots",showWarnings = FALSE, recursive = TRUE)
# plot power distributin vs model fit for scale-free topp. 
png(paste0("outputData/WGCNA/plots/power-vs-model-fit.png"), width=5*ppi, height=5*ppi, res=ppi)
# pdf(file = "outputData/WGCNA/plots/power-vs-model-fit.pdf", width = 5, height = 5, pointsize = 1)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")  
abline(h=0.9,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power. Again, we are looking at when the plot begins leveling out.
png(paste0("outputData/WGCNA/plots/median-connectivity.png"), width=5*ppi, height=5*ppi, res=ppi)
# pdf(file = "outputData/WGCNA/plots/median-connectivity.pdf", width = 5, height = 5, pointsize = 1)
plot(sft$fitIndices[,1], sft$fitIndices[,6],
     xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",
     main = paste("Median connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
dev.off()

####### set threshold power 
softThres <- 20

# generate adjacency matrix for specific soft threshopld 
ADJ1=abs(cor(norm_counts, use="p"))^softThres

# calculate connectivity 
k=softConnectivity(datE=norm_counts, power=softThres)

# Plot a histogram of k 
#png(paste0("outputData/plots/WGCNAplots/degree-dist-", softThres, ".png"), width=5*ppi, height=5*ppi, res=ppi)
pdf(file = paste0("outputData/WGCNA/plots/degree-dist-", softThres, ".pdf"), width = 5, height = 5, pointsize = 1)
hist(k, breaks = 50, col = "indianred", )
dev.off()

# plot scale free topo. fit for soft threshold power 
#png(paste0("outputData/plots/WGCNAplots/model-fit-soft-threshold-", softThres, ".png"), width=5*ppi, height=5*ppi, res=ppi)
pdf(file = paste0("outputData/WGCNA/plots/model-fit-soft-threshold-", softThres, ".pdf"), width = 5, height = 5, pointsize = 1)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run WGCNA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# assign cor function to WGCNA 
cor <- WGCNA::cor

# cut height 0.25 
net_0.25 <- blockwiseModules(norm_counts,
                          networkType = "signed",
                          maxBlockSize = 30000, 
                          TOMType = "unsigned", 
                          power = softThres, 
                          numericLabels = TRUE, 
                          randomSeed = 1234,
                          mergeCutHeight = 0.25,
                          minModuleSize = 30,
                          deepSplit = 2)

##### reassign core function to stats package 
cor<-stats::cor

# save network to file 
saveRDS(net_0.25, file = "outputData/WGCNA/network-signed_cutHeight-0.25.rds")

# save module eigenegenes 
saveRDS(net_0.25$MEs, file = paste0("outputData/WGCNA/module_eigengenes_cutHeight-0.25.rds"))

# check number of modules identified 
table(net_0.25$colors)


dir.create("./outputData/WGCNA/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(table(net_0.25$colors)), file = paste0("outputData/WGCNA/tables/gene-counts-per-module-network-0.25_colors.csv"))

# calculate kme for all nodes 
datKME=signedKME(norm_counts, net_0.25$MEs, outputColumnName="MM.")
saveRDS(datKME, file = paste0("outputData/WGCNA/kME-signed.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize network 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract genes that were assigned to modules 
restGenes = names(net_0.25$colors[!net_0.25$colors==0])

# colors for genes that were assigned to modules 
moduleLabels = net_0.25$colors
moduleColors = labels2colors(moduleLabels)
restGenes_cols = moduleColors[!net_0.25$colors==0]

# calc adjacency matrix for genes assigned to modules 
adjacency <- adjacency(norm_counts[, colnames(norm_counts) %in% restGenes], power = softThres , type = "signed") 

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM <- 1-TOM

# perform hierarchical clustering of the TOM
hier1 = hclust(as.dist(dissTOM), method="average" )

# set diagonal to NA
diag(dissTOM) = NA

# plot tree and module assigned from h.clustering of TOM 
ppi=1200
png(paste0("outputData/WGCNA/plots/dendrogram.png"), width=10*ppi, height=5*ppi, res = ppi, bg = "transparent")

plotDendroAndColors(hier1, restGenes_cols, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# plot TOM heatmap based on topological overlap 
ppi=1200

png(paste0("outputData/WGCNA/plots/TOM-plot.png"), width=5*ppi, height=5*ppi, res=ppi, bg = "transparent")
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(dissTOM^10, hier1, main = "TOM heatmap plot, module genes" , Colors = restGenes_cols, col=myheatcol)
dev.off()

####### visualize clustering of module eigengenes ########

# remove unassigned module 
ME_sub <- net_0.25$MEs[, !c(colnames(net_0.25$MEs) %in% "ME0")]

# calculate dissamility matrix 
dissimME=(1-cor(ME_sub))

# cluster using distances between dissimilarity of eigenegenes 
hclustdatME=hclust(as.dist(dissimME), method="average" )

# Plot the eigengene dendrogram

pdf(file = paste0("outputData/WGCNA/plots/dissimilarity-based-ME-clustering.pdf"), width=7, height=5, pointsize = 1)
plot(hclustdatME, main="Clustering tree based on module eigengene similarity", las = 1)
abline(h=0.25, col = "red") 
abline(h=0.3, col = "red") 
abline(h=0.4, col = "red") 
dev.off() 

# plot pairwise correlations 

pdf(file = paste0("outputData/WGCNA/plots/dissimilarity-based-ME-clustering-cors.pdf"), width=15, height=15, pointsize = 1)
plotMEpairs(net_0.25$MEs)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save sessionInfo 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save sessInfo 
writeLines(capture.output(sessionInfo()), "outputData/WGCNA/WGCNA_sessionInfo.txt")

