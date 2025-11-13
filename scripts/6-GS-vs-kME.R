
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Compare Gene significance scores vs module connectivity measure (kME)
#
# Lab: Cramer 
# Written by Owen Wilkins
# Edited by Charles Puerner
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(WGCNA)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggrepel)

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

# load network
net <- readRDS(paste0("outputData/WGCNA/network-signed_cutHeight-0.25.rds"))

dir.create("./outputData/WGCNA/GS_vs_kME_data", showWarnings = F, recursive = T)

# load module eigengenes 
kme <- readRDS(paste0("outputData/WGCNA/kME-signed.rds"))
write.csv(kme, file = "outputData/WGCNA/GS_vs_kME_data/signed-kME-matrix.csv")

# read in DE results for GS
deg <- read.csv("outputData/DEG/biofilm_vs_planktonic_limma-DEGs_all-genes.csv", stringsAsFactors = FALSE, row.names = 1)

# total adjacency 
total_adj <- read.csv(paste0("outputData/WGCNA/network/total-adjacency.csv"))

# degree df 
degree <- read.csv(paste0("outputData/WGCNA/network/network-degree-tom-thres-0.95-percentile.csv"))
degree <- degree[, c(2, 3)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize GS vs kme for modules of interest 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Create directories to deposit the data
  dir.create(paste0("outputData/WGCNA/GS_vs_kME_data/connectivity/"), showWarnings = TRUE, recursive = FALSE)
  dir.create(paste0("outputData/WGCNA/GS_vs_kME_data/connectivityPlots/"), showWarnings = TRUE, recursive = FALSE)
  
  # calculate cor matrix 
  ADJ1 <- abs(cor(norm_counts, use="p"))^20
  
  # confirm cors in same order as names of modules 
  identical(colnames(ADJ1), names(net$colors))
  
  # calc intramodular connectivity 
  int_conn <- intramodularConnectivity(ADJ1, net$colors)
  
  # remove ME0 genes 
  me0_genes <- names(which(net$colors=="0"))
  int_conn <- int_conn[-which(rownames(int_conn) %in% me0_genes),]
  kme <- kme[-which(rownames(kme) %in% me0_genes),]
  
  # filter out deg genes not included in kme df (module 0 genes)
  deg <- deg[rownames(deg) %in% rownames(kme),]
  
  # filter out low expressed genes not tested for DEG from kme 
  kme <- kme[rownames(kme) %in% rownames(deg),]
  int_conn <- int_conn[rownames(int_conn) %in% rownames(deg),]
  degree <- degree[degree$gene %in% rownames(deg),]
  total_adj <- total_adj[total_adj$gene %in% rownames(deg),]
  
  # add in genes with zero degree so that it can be combined with other dfs 
  missing_genes <- total_adj[!total_adj$gene %in% degree$gene,]
  missing_genes_df <- data.frame(gene = missing_genes$gene, degree = 0)
  degree <- rbind(degree, missing_genes_df)
  
  # order kme, int_conn, & deg same way 
  deg <- deg[order(rownames(deg)),]
  kme <- kme[order(rownames(kme)),]
  degree <- degree[order(degree$gene),]
  total_adj <- total_adj[order(total_adj$gene),]
  int_conn <- int_conn[order(rownames(int_conn)),]
  
  # confirm identical ordering 
  identical(rownames(deg), rownames(kme))
  identical(rownames(deg), rownames(int_conn))
  identical(rownames(deg), degree$gene)
  identical(rownames(deg), total_adj$gene)
  
  # function to plot GS vs kme and export info for modules of interest 
  plot_GS_vs_kme <- function(module){
    message("running module", module)
    # generate df of data of interest 
    df_tmp_GS_kme <- data.frame("gene" = rownames(deg), "GS" = -log10(deg$adj.P.Val), "kme" = kme[, c(paste0("MM.", module))],
                                "kTotal" = int_conn$kTotal, "kWithin" = int_conn$kWithin, "kOut" = int_conn$kOut, "kDiff" = int_conn$kDiff, 
                                "DE-padj" = deg$adj.P.Val, "log2FoldChange" = deg$logFC, "degree_at_TOM_95p" = degree$degree, "total_adjacency_within_module" = total_adj$total_adjacency)
    # get genes in module 
    module_genes <- names(net$colors)[net$colors == module]
    # subset to module genes only 
    df_tmp_GS_kme_sub <- df_tmp_GS_kme[df_tmp_GS_kme$gene %in% module_genes,]
    # order by kme 
    df_tmp_GS_kme_sub <- df_tmp_GS_kme_sub[order(df_tmp_GS_kme_sub$kme, decreasing=TRUE),]
    # write to file  
    write.csv(df_tmp_GS_kme_sub, file = paste0("./outputData/WGCNA/GS_vs_kME_data/connectivity/connectivity-characteristics-module-", module, ".csv"), row.names = FALSE)
    
    # plot GS vs kme 
   
    pdf(paste0("./outputData/WGCNA/GS_vs_kME_data/connectivityPlots/GS-vs-kme-module-", module, ".pdf"), width=7, height=6.5)
    p <- ggscatter(df_tmp_GS_kme_sub, x = "kme", y = "GS",
              color = "#D3352B", shape = 1, size = 2.5,
              xlab = "kME (Module membership)",
              ylab = "Gene significance (-log10 P-value)") +
      ggtitle(paste0("Module ", module)) + 
      geom_label_repel(data = subset(df_tmp_GS_kme_sub, GS > max(df_tmp_GS_kme_sub$GS)*0.9 & kme > 0.8), 
                       aes(label = gene), 
                       point.padding = 0.5,
                       label.size = 0.1,
                       segment.size = 0.3,
                       segment.color = 'grey50', size = 3) +
      theme(plot.title = element_text(size = 20))
    print(p)
    dev.off()
  }

  # apply function to modules of interest 
  module_list <- sapply(colnames(ME_sub), function(x) strsplit(x, "E")[[1]][2])
  lapply(module_list, plot_GS_vs_kme)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save sessionInfo 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save sessInfo 
writeLines(capture.output(sessionInfo()), "./outputData/WGCNA/GS_vs_kME_data/GS_vs_kME_sessionInfo.txt")
