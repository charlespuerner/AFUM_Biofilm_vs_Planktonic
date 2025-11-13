#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Title: Mixed effects DEG 

# Description: Mixed effects differential expression analysis of biofilm vs planktonic growth in normoxia conditions
#
# Lab: Cramer 
# Author: Owen Wilkins
# Edited by: Charles Puerner
#
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
library(rtracklayer)
library(edgeR)
library(statmod)
library(circlize)


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
identical(colnames(counts), rownames(meta))

# generate DGE oject containing counts and sample data 
dge <- DGEList(counts= counts, samples = meta)

# check dist of lib sizes 
hist(dge$samples$lib.size/1e6)

# ratio from smallest lib size to largest
max(dge$samples$lib.size)/min(dge$samples$lib.size)

# set model matrix 
design = model.matrix( ~culture, meta)

#filter out genes with consistently low or zero counts
keep <- filterByExpr(dge, design = design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#normalize the counts using the trimmed mean of M-vlaues (TMM)
normdata <- calcNormFactors(dge, method = "TMM")

#first voom normalization
v <- voom(normdata, design, plot=TRUE)

# delineate that differences in samples are due to random effects
corfit <- duplicateCorrelation(v, design, block=meta$condition)
corfit$consensus.correlation

#run voom a second time with the random effects accounted for
vobj = voom(normdata, 
            design, 
            plot=TRUE, 
            block=meta$condition, 
            correlation=corfit$consensus)

# Estimate linear mixed model with a single variance component
corfit <- duplicateCorrelation(vobj, design, block=meta$condition)
corfit$consensus.correlation

# But this step uses only the genome-wide average for the random effect
fit <- lmFit(vobj, design, block=meta$condition, correlation=corfit$consensus)

# Fit Empirical Bayes for moderated t-statistics
fit <- eBayes(fit, robust = FALSE, trend = FALSE)

# check test results available 
summary(decideTests(fit))

# return results table 
res <- topTable(fit, coef = "culturebiofilm", n = Inf, sort.by = "p")
res <- res[order(res$adj.P.Val),]

# check how many genes sign. at thresholds 
table(res$adj.P.Val<0.05)
table(res$adj.P.Val<0.05 & res$logFC<0)
table(res$adj.P.Val<0.05 & res$logFC < -2)
table(res$adj.P.Val<0.05 & abs(res$logFC)> 2)

# quick volcano plot
plot(res$logFC, -log10(res$adj.P.Val))
abline(h=-log10(0.05))

# write results to file 
res_sub <- res[res$adj.P.Val<0.05,]
write.csv(res_sub, paste0("outputData/DEG/biofilm_vs_planktonic_limma-DEGs_FDR-0.05.csv"))
write.csv(res, paste0("outputData/DEG/biofilm_vs_planktonic_limma-DEGs_all-genes.csv"))

# extract normalized counts 
norm_counts <- cpm(dge, log=FALSE, prior.count=3)
write.csv(norm_counts, file = "./outputData/DEG/CPM_biofilm_vs_planktonic.csv" )
norm_counts <- cpm(dge, log=TRUE, prior.count=3)
write.csv(norm_counts, file = "./outputData/DEG/CPM_biofilm_vs_planktonic_log.csv" )
# plot genes of interest 
boxplot(norm_counts['Afu4g06290',] ~ dge$samples$culture)
res_sub['Afu4g06290',]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# volcano plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load custom volcano plotting function 
source(paste0("./scripts/volcano.R"))

#add gene to results as column 
res$gene <- rownames(res)

# make volcano plot 
volcano_plot(dat = res, 
             out_dir = paste0("outputdata/DEG/plots/"), 
             out_name = paste0("volcano.pdf"), 
             title = paste0(" Differential gene expression"),
             alpha=0.05, fc_line = 1.5,
             alpha_label = 1e-7, 
             fc_label=2.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# heatmap 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select top X no. of variable genes 
DEgenes_df <- res_sub[abs(res_sub$logFC)>2.5 & res_sub$adj.P.Val<0.05,]
DEgenes <- rownames(DEgenes_df)

# set up gene expression matrix 
mat1 <- norm_counts[DEgenes,]

# scale matrix by each col. values 
mat_scaled = t(apply(mat1, 1, scale))
colnames(mat_scaled) <- colnames(mat1)

# set up colors for heatmap 
cols = colorRamp2(c(-2, 0, 2), c("#4575B4", "#FFFFBF", "#D73027"))
cols1 <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(12, "Paired")
cols3 <- brewer.pal(11, "RdYlBu")


# set up annotation bar for samples 
ha1 = HeatmapAnnotation(Culture_condition = as.factor(dge$samples$culture), 
                        col = list(Culture_condition = c("biofilm" = cols3[3], "planktonic" = cols3[9]), show_legend = TRUE))

ha2 = HeatmapAnnotation(Growth = as.factor(dge$samples$growth), 
                        col = list(Growth = c("biofilm" = cols1[2], 
                                               "normoxia_planktonic" = cols1[3],
                                               "hypoxia_planktonic" = cols1[1]), 
                                   show_legend = TRUE))

### strain anno 
ha3 = HeatmapAnnotation(Strain = as.factor(dge$samples$strain), 
                        col = list(Strain = c("AF293" = cols1[4],"CEA10" = cols1[5]), show_legend = TRUE))

ha4 = HeatmapAnnotation(Condition = as.factor(dge$samples$condition), 
                        col = list(Condition = c("hypoxia" = cols2[3], "normoxia" = cols2[7]), show_legend = TRUE))

ha = columnAnnotation(x = anno_text(rownames(dge$samples), 
                                    which="column", 
                                    rot = 45, 
                                    gp = gpar(fontsize = 12)))


# generate heatmap object Mi
ht1 = Heatmap(mat_scaled, name = "Expression", #col = cols, 
              top_annotation = c(ha1, ha2, ha3, ha4),
              bottom_annotation = ha,
              show_row_names = FALSE, show_column_names = FALSE)

# plot the heatmap 

pdf(paste0("outputData/DEG/plots/heatmap_DEG_limma_padj-0.05_lfc-2.5_biofilm", ".pdf"), width=10, height=8, pointsize = 1, bg = "transparent")
draw(ht1, row_title = "Genes", column_title = "DEGs - P < 0.05, abs. LFC > 2.5")
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MA plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(paste0("scripts/ma-plot.R"))

# make ma plot 
ma_plot(dat = res, 
               out_dir = paste0("outputdata/DEG/plots/"), 
               out_name = paste0("ma.pdf"), 
               title = paste0(" Differential gene expression"),
               alpha=0.05, 
               alpha_label = 1e-40, 
               ylimit_low=-5, ylimit_up = 5,
               fc_label=3, mean_label = 1)
