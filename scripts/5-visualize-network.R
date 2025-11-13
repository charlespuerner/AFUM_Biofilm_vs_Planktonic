
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: Network visualization & hub gene identification

# Description: Perform basic network vis. of network connections, and hub gene idenification analyses 
# 			   (hub gene ID analyses: TOM threshold downsampling sensitivity analysis, adjacency-based) 
#
# Lab: Cramer 
# Author: Owen Wilkins
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
library(ggrepel)
library(RCy3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setwd 
setwd("/Users/f0052zp/Documents/Research_data_DartmouthDrive/Genomics/Bioinformatic_Rdata/Analysis_Normoxia_Biofilm/biofilm_vs_norm-and-lowO2_planktonic")

# read in normalized counts 
norm_counts <- readRDS("outputData/Exploratory/normalized_counts_vst.rds")

# read in phenotypic data 
meta <- read.csv("inputData/metadata.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
# meta$sample <- NULL
meta$growth <- factor(meta$growth, levels = c("normoxia_planktonic","hypoxia_planktonic","biofilm"))
meta$condition <- factor(meta$condition, levels = c("normoxia","hypoxia"))
meta$batch_no <- factor(as.character(meta$batch_no), levels = c("7",  "2",  "3",  "4",  "1",  "B1", "B2", "B3", "B4"))
meta$group <- factor(meta$group, levels = unique(meta$group))
meta$culture <- factor(meta$culture, levels = (c("planktonic","biofilm")))

# load module eigengenes 
MEs <- readRDS(paste0("outputData/WGCNA/module_eigengenes_cutHeight-0.25.rds"))

# remove unassigned module 
ME_sub <- MEs[, !c(colnames(MEs) %in% "ME0")]

# load module eigengenes 
net <- readRDS(paste0("outputData/WGCNA/network-signed_cutHeight-0.25.rds"))

# load kme
kme <- readRDS(paste0("outputData/WGCNA/kME-signed.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize network connections per module 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.create("./outputData/WGCNA/network/network_plots", showWarnings = F, recursive = T)

module_list <- sapply(colnames(ME_sub), function(x) strsplit(x, "E")[[1]][2])

# load custom functions for network analysis 
source("./Scripts/WGCNA/network-analysis-utilities.R")

# apply to modules of interest 
networks <- lapply(module_list, build_network, "#E98113", 100, 0.95)

# spot checks (compare visually to graphs and confirm correct values being plotted at .95 tom thres)
degree(networks$ME9)[order(degree(networks$ME9))]
degree(networks$ME3)[order(degree(networks$ME3))]
degree(networks$ME2)[order(degree(networks$ME2))]

# create df of degree for use in output spreadsheet
net_degree <- lapply(networks, function(x) degree(x)[order(degree(x), decreasing = TRUE)])


# organize list into single df 
degree_df_list <- list()
for(i in 1:length(net_degree)){
  degree_df_list[[i]] <- data.frame(gene = names(net_degree[[i]]), degree = net_degree[[i]])
  rownames(degree_df_list[[i]]) <- NULL
}
degree_df <- do.call(rbind, degree_df_list)
rownames(degree_df) <- NULL

# write to file 
write.csv(degree_df, file = paste0("./outputData/WGCNA/network/network-degree-tom-thres-0.95-percentile.csv"))

# plot the top genes from the networks 
mapply(plot_network, networks, names(networks), 50, 100) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# perform TOM thresholding sensitivity analysis 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("./scripts/WGCNA/sensitivity-analysis-function.R")

# run sensitivity analysis across various TOM thresholds 
top_genes_list <- lapply(module_list, sensitivity_analysis_degree, 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# adjacency-based hub gene identification
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("./scripts/WGCNA/total-adjacency-analysis-function.R")

dir.create("./outputData/WGCNA/network/networkSummarySheets", showWarnings = F, recursive = T)

# run total adjacency analysis within each module 
total_adjacency <- list()
for(i in 1:length(module_list)){
  total_adjacency[[i]] <- total_adjacency_degree_analysis(module_list[[i]], 20, top_genes_list[[i]])
}

# rename 
names(total_adjacency) <- names(top_genes_list)

# output df for spreadsheet 
total_adjacency_df <- do.call(rbind, total_adjacency)
total_adjacency_df <- total_adjacency_df[, c(1,2)]
rownames(total_adjacency_df) <- NULL
write.csv(total_adjacency_df, file = paste0("outputData/WGCNA/network/total-adjacency.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TOM heatmaps 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# apply to modules of interest 
lapply(module_list, module_heatmap, 100, 20)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save sessionInfo 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save sessInfo 
writeLines(capture.output(sessionInfo()), "./outputData/WGCNA/network/sessionInfo.txt")

