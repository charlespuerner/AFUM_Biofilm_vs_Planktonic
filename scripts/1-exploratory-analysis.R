#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Process raw counts and perform exploratory data analysis 
#
# Lab: Cramer 
#
# Original Author: Owen Wilkins
#
# Edited by: Charlie Puerner
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(WGCNA)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(DESeq2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#setwd in working directory.

# read in counts matrix 
counts <- read.table("./inputData/biofilm_vs_planktonic.txt", sep = "\t", header=TRUE, stringsAsFactors = F)

# convert gene ID to rownames 
rownames(counts) <- counts$gene_ID
counts$gene_ID <- NULL

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
counts <- counts[, colnames(counts) %in% rownames(meta)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prep. & normalize dataset 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check order of samples 
identical(colnames(counts), meta$sample)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = meta, 
  design = ~1 
)

# run the DEseq2 analysis 
dds <- DESeq(dds)

# check size factors 
barplot(sizeFactors(dds), las = 1, col = 'indianred')

# drop genes with low counts 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# normalize 
vst <- vst(dds)
colData(vst) <- colData(dds)

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(vst) %>%
  t() # Transpose this data

saveRDS(normalized_counts, file = "outputData/Exploratory/normalized_counts_vst.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample distance and dimensionality analysis using pheatmap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### sample distance (euclidean) plot ###### 
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rownames(colData(vst)), sep = " - " )
colnames(sampleDistMatrix) <- paste( rownames(colData(vst)), sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)

##### 
# set up colors for heatmap 
cols = colorRampPalette( rev(brewer.pal(9, "RdPu")) )(255)
cols1 <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(12, "Paired")
cols3 <- brewer.pal(11, "RdYlBu")
# set up annotation bar for samples 
ha1 = HeatmapAnnotation(growth = as.factor(colData(dds)$growth),
                        col = list(growth = c("biofilm" = cols1[2], 
                                              "normoxia_planktonic" = cols1[3],
                                              "hypoxia_planktonic" = cols1[1]), show_legend = TRUE))

ha1.5 = HeatmapAnnotation(Culture_condition = as.factor(colData(dds)$culture),
                        col = list(Culture_condition = c("biofilm" = cols3[3], "planktonic" = cols3[9]), show_legend = TRUE))

# set up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData(dds)$sample,
                                   which="column",
                                   rot = 45,
                                   gp = gpar(fontsize = 10)))

### strain anno 
ha2 = HeatmapAnnotation(Strain = as.factor(meta$strain), 
                        col = list(Strain = c("AF293" = cols1[4],"CEA10" = cols1[5]), show_legend = TRUE))

ha3 = HeatmapAnnotation(Condition = as.factor(colData(dds)$condition), 
                                          col = list(Condition = c("hypoxia" = cols2[3], "normoxia" = cols2[7]), show_legend = TRUE))

# generate heatmap object Mi
ht1 = Heatmap(sampleDistMatrix, name = "sampleDistMatrix", col = cols, 
              top_annotation = c(ha3, ha1.5, ha1, ha2),
              height = 5,
              width = 2,
              bottom_annotation = ha,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 10),
              clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
              show_column_names = FALSE, )

# plot the heatmap 
pdf(paste0("outputData/Exploratory/plots/heatmap_euclidean_distance_matrix.pdf"), width=10, height=8, pointsize = 1, bg = "transparent")
draw(ht1, column_title = "Euclidean distance based hierarchical clustering",auto_adjust = F, background = "transparent", padding = unit(c(0, 0.25, 0, 0), "in"))
dev.off()

###### principal components analysis ###### 

# calculate gene expression level variance between samples (to determine how 
# many genes to use for dimensionality reduction and clustering)
var <- rev(rowVars(assay(vst))[order(rowVars(assay(vst)))])

# plot variance for genes across samples

pdf(paste0("outputData/Exploratory/plots/per_gene_variance.pdf"), width=5, height=5, pointsize = 1, bg = "transparent")
plot(var, las = 1, main="Sample gene expression variance", xlab = "Gene", ylab = "Variance")
abline(v=4000, col="red")

dev.off()

# perform PCA and order by variance 
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(4000, length(rv)))]
pca <- prcomp(t(assay(vst)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
names(percentVar)[1:10] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC19", "PC10")

# plot variance for top 10 PCs 
pdf(paste0("outputData/Exploratory/plots/PCA_variance.pdf"), width=5, height=5, pointsize = 1, bg = "transparent")
barplot(percentVar[1:10], col = "indianred", las = 1, ylab = "Percent Variance explained", cex.lab = 1.2, ylim =c(0,0.5))
dev.off()

# group into data frame and add sample labels 
pca_df <- as.data.frame(pca$x)
pca_df$growth <- as.factor(dds$growth)
pca_df$strain <- as.factor(dds$strain)
pca_df$sample_ids <- colnames(dds)
pca_df$batch <- dds$batch_no
pca_df$group <- dds$group

# add colors for plotting to df 
pca_df$col <- NA
for(i in 1:length(levels(pca_df$growth))){
  ind1 <- which(pca_df$growth == levels(pca_df$growth)[i])
  pca_df$col[ind1] <- i
}

pca_df$shape <- NA
for(i in 1:length(levels(pca_df$strain))){
  ind1 <- which(pca_df$strain == levels(pca_df$strain)[i])
  shapes <- c(16, 17) [i]
  pca_df$shape[ind1] <- shapes
}

pca_df$col_color <- NA
for (i in 1:length(levels(pca_df$growth))) {
  ind <- which(pca_df$growth == levels(pca_df$growth)[i])
  color_value <- c("#7570B3","#1B9E77","#D95F02")[i]
  pca_df$col_color[ind] <- color_value
}


# PC1 vs PC2
ppi = 300
# pdf(paste0("outputData/plots/pca_PC1vsPC2.png"), width=7, height=7, pointsize = 1, bg = "transparent")
png(paste0("outputData/Exploratory/plots/pca_PC1vsPC2.png"), width=7*ppi, height=7*ppi, res=ppi)
plot(pca_df[, 1], pca_df[, 2], 
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=(paste0("PC1 vs PC2 for", " most variable genes")),
     pch=pca_df$shape, 
     cex=2, cex.lab=1.3, cex.axis = 1.15, las=1, panel.first = grid(),
     col=pca_df$col_color,
     xlim=c(-55,100),
     ylim=c(-50, 50))

legend(15, 50, c(levels(pca_df$growth)), pch = 16, col = c("#7570B3","#1B9E77","#D95F02"), cex = 1)
legend(75 , 50,  c(levels(pca_df$strain)), pch = unique(pca_df$shape), cex = 1)
dev.off()



# PC1 vs PC3
png(paste0("outputData/Exploratory/plots/pca_PC1vsPC3.png"), width=7*ppi, height=7*ppi, res=ppi)
plot(pca_df[, 1], pca_df[, 3], 
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
     ylab = paste0("PC3 (", (round(percentVar[3], digits=3)*100), "% variance)"),
     main=(paste0("PC1 vs PC3 for", " most variable genes")),
     pch=pca_df$shape, 
     cex=2, cex.lab=1.3, cex.axis = 1.15, las=1, panel.first = grid(),
     col=pca_df$col_color,
     xlim=c(-60, 100),
     ylim=c(-40, 65))

legend(15, 60, c(levels(pca_df$growth)), pch = 16, col = c("#7570B3","#1B9E77","#D95F02"), cex = 1)
legend(75 , 60,  c(levels(pca_df$strain)), pch = unique(pca_df$shape), cex = 1)
dev.off()

# PC2 vs PC3
png(paste0("outputData/Exploratory/plots/pca_PC2vsPC3.png"), width=7*ppi, height=7*ppi, res=ppi)
plot(pca_df[, 3], pca_df[, 2], 
     xlab = paste0("PC3 (", (round(percentVar[3], digits=3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=(paste0("PC3 vs PC2 for", " most variable genes")),
     pch=pca_df$shape, 
     cex=2, cex.lab=1.3, cex.axis = 1.15, las=1, panel.first = grid(),
     col=pca_df$col_color,
     xlim=c(-40, 65),
     ylim=c(-50, 50))

legend(10, -36, c(levels(pca_df$growth)), pch = 16, col = c("#7570B3","#1B9E77","#D95F02"), cex = 1)
legend(50 , -36,  c(levels(pca_df$strain)), pch = unique(pca_df$shape), cex = 1)
dev.off()


library(rgl)
plot3d(pca_df[, 1], pca_df[, 2],pca_df[, 3],
       col=pca_df$col_color,type = "s", size = 3,
       xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
       ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
       zlab = paste0("PC3 (", (round(percentVar[3], digits=3)*100), "% variance)"),
       rgl.printRglwidget = TRUE)
grid3d(c("x", "y+", "z"))
rgl.postscript("./outputData/Exploratory/plots/3Dpca_PC1_PC2_PC3.pdf", fmt = 'pdf')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hierachical clustering analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select top X no. of variable genes 
topVarGenes <- head(order(rowVars(assay(vst)), decreasing=TRUE), 4000)

# set up gene expression matrix 
mat1 <- assay(vst)[topVarGenes,]

# scale matrix by each col. values 
mat_scaled = t(apply(mat1, 1, scale))
colnames(mat_scaled) <- colData(dds)$sample

# set up colors for heatmap 
cols = colorRamp2(c(-3,0, 3), c("#4575B4", "#FFFFBF", "#D73027"))
cols1 <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(12, "Paired")
cols3 <- brewer.pal(11, "RdYlBu")
# set up annotation bar for samples 
ha1 = HeatmapAnnotation(growth = as.factor(colData(dds)$growth),
                        col = list(growth = c("biofilm" = cols1[2], 
                                              "normoxia_planktonic" = cols1[3],
                                              "hypoxia_planktonic" = cols1[1]), show_legend = TRUE))

ha1.5 = HeatmapAnnotation(Culture_condition = as.factor(colData(dds)$culture),
                          col = list(Culture_condition = c("biofilm" = cols3[3], "planktonic" = cols3[9]), show_legend = TRUE))

# set up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData(dds)$sample,
                                    which="column",
                                    rot = 45,
                                    gp = gpar(fontsize = 4)))

### strain anno 
ha2 = HeatmapAnnotation(Strain = as.factor(meta$strain), 
                        col = list(Strain = c("AF293" = cols1[4],"CEA10" = cols1[5]), show_legend = TRUE))

ha3 = HeatmapAnnotation(Condition = as.factor(colData(dds)$condition), 
                        col = list(Condition = c("hypoxia" = cols2[3], "normoxia" = cols2[7]), show_legend = TRUE))

# generate heatmap object Mi
ht1 = Heatmap(mat_scaled, name = "Expression", #col = col, 
              top_annotation = c(ha3, ha1.5, ha1, ha2),
              #bottom_annotation = ha,
              show_row_names = FALSE, show_column_names = FALSE)

# plot the heatmap 

pdf(paste0("outputData/Exploratory/plots/heatmap_most-variable-4k-genes.pdf"), width=10.5, height=6, pointsize = 1, bg = "transparent")
draw(ht1, row_title = "Genes", column_title = "Unsupervised hierarchical clustering")
dev.off()


