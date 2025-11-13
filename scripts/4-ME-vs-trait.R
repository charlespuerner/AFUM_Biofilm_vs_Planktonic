
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: Module vs trait (normoxic/hypoxic)

# Description: Mixed effects analysis to identify module-trait associations for normoxia/hypoxia tx condition
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
library(limma)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setwd 
setwd("/Users/f0052zp/Documents/Research_data_DartmouthDrive/Genomics/Bioinformatic_Rdata/Analysis_Normoxia_Biofilm/biofilm_vs_norm-and-lowO2_planktonic")

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
meta$strain <-factor(meta$strain, levels = unique(meta$strain))

# remove samples failing QC # remove samples failing QC 
meta <- meta[!meta$qc == "FAIL",]

# subset counts for samples passing QC 
norm_counts <- norm_counts[rownames(norm_counts) %in% rownames(meta),]

# load module eigengenes 
MEs_0.25 <- readRDS(paste0("outputData/WGCNA/module_eigengenes_cutHeight-0.25.rds"))

# remove unassigned module 
MEs_0.25_sub <- MEs_0.25[, !c(colnames(MEs_0.25) %in% "ME0")]

# load network 
net_0.25 <- readRDS(paste0("outputData/WGCNA/network-signed_cutHeight-0.25.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test Module vs trait relations 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define numbers of genes and samples
nGenes = ncol(norm_counts)
nSamples = nrow(norm_counts)


# confirm samples are in same order in meta and ME df 
identical(rownames(meta), rownames(MEs_0.25_sub))

# transpose modules into rows  
MEs_0.25_sub_t <- t(MEs_0.25_sub)

# set model matrix 
design = model.matrix( ~ 0 + culture, meta)

# delineate that differences in patient samples are due to random effects
corfit <- duplicateCorrelation(MEs_0.25_sub_t, design, block=meta$strain)
corfit$consensus.correlation

# But this step uses only the genome-wide average for the random effect
fit <- lmFit(MEs_0.25_sub_t, design, block=meta$strain, correlation=corfit$consensus)

# Fit Empirical Bayes for moderated t-statistics
fit <- eBayes(fit, robust = FALSE, trend = FALSE)

# check test results available 
summary(decideTests(fit))

# return results tadesign# return results table Clade 1


res <- topTable(fit, coef = "culturebiofilm", n = Inf, sort.by = "p")
res <- res[order(res$adj.P.Val),]

dir.create("./outputData/WGCNA/ME_vs_trait")

# write results to file 
res_sub <- res[res$adj.P.Val<0.05,]
write.csv(res_sub, paste0("outputData/WGCNA/ME_vs_trait/limma-differential-module-activation-", "culturebiofilm", "FDR-0.05.csv"))
write.csv(res, paste0("outputData/WGCNA/ME_vs_trait/limma-differential-module-activation-", "culturebiofilm", "all-modules.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize Module vs trait relations 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####### heatmap of samples vs eigenvalues/eigengenes #######

dir.create("./outputData/WGCNA/ME_vs_trait/plots", showWarnings = F, recursive = T)

# add sample names as column to meta 
meta$sample <- rownames(meta)

# set up colors for heatmap 
cols = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cols1 <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(12, "Paired")
cols3 <- brewer.pal(11, "RdYlBu")

# set up annotation bar for samples 

### strain anno 
ha1 = HeatmapAnnotation(Culture_condition = as.factor(meta$culture), 
                        col = list(Culture_condition = c("biofilm" = cols3[3], "planktonic" = cols3[9]), show_legend = TRUE))

ha2 = HeatmapAnnotation(Growth = as.factor(meta$growth), 
                        col = list(Growth = c("biofilm" = cols1[2], 
                                              "normoxia_planktonic" = cols1[3],
                                              "hypoxia_planktonic" = cols1[1]), 
                                   show_legend = TRUE))

ha3 = HeatmapAnnotation(Strain = as.factor(meta$strain), 
                        col = list(Strain = c("AF293" = cols1[4],"CEA10" = cols1[5]), show_legend = TRUE))

ha4 = HeatmapAnnotation(O2_Condition = as.factor(meta$condition),
                        col = list(O2_Condition = c("hypoxia" = cols2[3], "normoxia" = cols2[7]), show_legend = TRUE))

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(meta$sample, which="column", rot = 45, gp = gpar(fontsize = 4)))

# generate heatmap for 0.25 cutheight network
ht1 = Heatmap(t(MEs_0.25_sub), name = "MEs", #col = col, 
              top_annotation = c(ha1, ha2, ha3),
              bottom_annotation = ha,
              show_row_names = T, show_column_names = F, row_names_gp = gpar(fontsize = 9))

# plot the heatmap 

pdf(paste0("outputData/WGCNA/ME_vs_trait/plots/heatmap_samples-vs-eigengenes-cutHeight-0.25.rds.pdf"), width=10.5, height=6, pointsize = 1)
draw(ht1, row_title = "Eigengenes", column_title = "Eigengene-based clustering (cutheight - 0.25)")
dev.off()


#######  heatmap & boxplots of specific eigengenes ####### 

# vectorize gene - module relations 
gene_module_key <- tibble::enframe(net_0.25$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))

# spot checks 
net_0.25$colors["Afu1g03992"]
gene_module_key[gene_module_key$gene=="Afu1g03992",]

dir.create("./outputData/WGCNA/ME_vs_trait/ME_heatmaps", showWarnings = F, recursive = T)
dir.create("./outputData/WGCNA/ME_vs_trait/ME_vs_trait_boxplots", showWarnings = F, recursive = T)
# function for module heatmap generation 
source("Scripts/WGCNA/module-heatmap-function.R")
source("Scripts/WGCNA/module-boxplot-function.R")

# generate heatmap for modules of interest 
lapply(rownames(MEs_0.25_sub_t), make_module_heatmap, module_eigengenes_df = MEs_0.25_sub)



# generate boxplot for modules of interest 
lapply(rownames(MEs_0.25_sub_t), make_module_boxplot, "culture", MEs_0.25_sub)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# strain specific trait association analysis 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vector of unique strains to test 
strains <- unique(meta$strain)

# empty matrices to fill with pvals and slopes 
res_mat_p <- matrix(NA, nrow = nrow(MEs_0.25_sub_t), ncol = length(strains))
res_mat_b <- matrix(NA, nrow = nrow(MEs_0.25_sub_t), ncol = length(strains))

# loop over strains and MEs and test for assocition with hypoxia status 

for(i in 1:length(strains)){
  # get sample ids of strain of interest
  strains_sub <- colnames(MEs_0.25_sub_t)[colnames(MEs_0.25_sub_t) %in% rownames(meta)[meta$strain==strains[i]]]
  # subset ME data for these samples
  MEs_0.25_sub_t_sub <- MEs_0.25_sub_t[, c(colnames(MEs_0.25_sub_t) %in% strains_sub)]
  # subset metadata for strain
  meta_sub <- meta[c(rownames(meta) %in% strains_sub), ]
  # if names of samples match, loop over MEs and test each for association within strain
  if(identical(rownames(meta_sub), colnames(MEs_0.25_sub_t_sub))){
    for(j in 1:nrow(MEs_0.25_sub_t_sub)){
      lm_tmp <- summary(lm(MEs_0.25_sub_t_sub[j,] ~ as.factor(meta_sub$culture)))
      res_mat_b[j, i] <- lm_tmp$coefficients[2,1]
      res_mat_p[j, i] <- lm_tmp$coefficients[2,4]
      rm(lm_tmp)
    }
  } else{message("FAIL, check input")}
  # adjust each column for multiple testing 
  res_mat_p[, i] <- p.adjust(res_mat_p[,i], method = "BH")
  # log transform
  res_mat_p[, i] <- -log10(res_mat_p[, i])
}

# add row and column names 
rownames(res_mat_p) <- rownames(MEs_0.25_sub_t)
rownames(res_mat_b) <- rownames(MEs_0.25_sub_t)
colnames(res_mat_p) <- strains
colnames(res_mat_b) <- strains

# set up colors for heatmap 
cols = colorRamp2(c(0, 5.5), c("#F5F3F3", "red"))
cols1 <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(12, "Paired")
cols3 <- brewer.pal(11, "RdYlBu")

# set up annotation bar for samples 

### strain anno 

ha = columnAnnotation(x = anno_text(colnames(res_mat_p), which="column", rot = 45, gp = gpar(fontsize = 9)))


# generate heatmap for 0.25 cutheight network
ht1 = Heatmap(res_mat_p, 
              name = "Significance", 
              col = cols, 
              bottom_annotation = ha ,
              show_row_names = T, 
              show_column_names = F, 
              row_names_gp = gpar(fontsize = 9), 
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(res_mat_p[i, j] > -log10(0.05)) {
                  grid.text("*", x, y)
                } 
              })


pdf(paste0("./outputData/WGCNA/ME_vs_trait/heatmap_ME-significnace.pdf"), width=4, height=8, pointsize = 1)
draw(ht1, row_title = "Eigengenes", column_title = "Strain-specific ME-biofilm association significance")
dev.off()


##### plot for betas 

# change color ramp for beta scale 
cols = colorRamp2(c(-0.5, 0, 0.5), c("blue", "#F5F3F3", "red"))

# generate beta heatmap 
ht1 <- Heatmap(res_mat_b, 
        name = "Direction", 
        col = cols, 
        bottom_annotation = ha,
        show_row_names = T, 
        show_column_names = F, 
        row_names_gp = gpar(fontsize = 9))

# 

pdf(paste0("outputData/WGCNA/ME_vs_trait/heatmap_ME-direction.pdf"), width=4, height=8, pointsize = 1)
draw(ht1, row_title = "Eigengenes", column_title = "Strain-specific ME-Biofilm association directions")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save sessionInfo 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save sessInfo 
writeLines(capture.output(sessionInfo()), "outputData/WGCNA/ME_vs_trait/ME_vs_trait_sessionInfo.txt")

