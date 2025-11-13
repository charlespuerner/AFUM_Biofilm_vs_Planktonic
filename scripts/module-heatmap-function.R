make_module_heatmap <- function(module_name,
                                expression_mat = norm_counts,
                                metadata_df = meta,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df) {
  
  # Set up the module eigengene with its sample
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample")
  
  # merge eigengene values into meta df 
  meta_merge <- merge(meta, module_eigengene, by = "sample")
  
  
  # get gene names of those in module 
  genes_to_keep <- gene_module_key$gene[gene_module_key$module == module_name]
  
  # subset counts for module genes 
  norm_counts_sub <- norm_counts[, colnames(norm_counts) %in% genes_to_keep]
  
  # scale matrix by each col. values 
  mat_scaled = t(apply(norm_counts_sub, 2, scale))
  colnames(mat_scaled) <- rownames(norm_counts_sub)
  
  mat_scaled <- mat_scaled[, order(colnames(mat_scaled))]
  meta_merge <- meta_merge[order(meta_merge$sample), ]
  
  #identical(colnames(mat_scaled), meta_merge$sample)
  
  
  cols = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  cols1 <- brewer.pal(8, "Dark2")
  cols2 <- brewer.pal(12, "Paired")
  cols3 <- brewer.pal(11, "RdYlBu")
  
  # set up annotation bar for samples 
  
  ### strain anno 
  ha1 = HeatmapAnnotation(Culture_condition = as.factor(meta_merge$culture), 
                          col = list(Culture_condition = c("biofilm" = cols3[3], "planktonic" = cols3[9]), show_legend = TRUE))
  
  ha1.5 = HeatmapAnnotation(module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(meta_merge, module_name)))
  
  ha2 = HeatmapAnnotation(Growth = as.factor(meta_merge$growth), 
                          col = list(Growth = c("biofilm" = cols1[2], 
                                                "normoxia_planktonic" = cols1[3],
                                                "hypoxia_planktonic" = cols1[1]), 
                                     show_legend = TRUE))
  
  ha3 = HeatmapAnnotation(Strain = as.factor(meta_merge$strain), 
                          col = list(Strain = c("AF293" = cols1[4],"CEA10" = cols1[5]), show_legend = TRUE))
  
  ha4 = HeatmapAnnotation(O2_Condition = as.factor(meta_merge$condition),
                          col = list(O2_Condition = c("hypoxia" = cols2[3], "normoxia" = cols2[7]), show_legend = TRUE))
  
  # se up column annotation labels (samples)
  ha = columnAnnotation(x = anno_text(meta_merge$sample, which="column", rot = 45, gp = gpar(fontsize = 4)))
  # generate heatmap object Mi
  ht1 = Heatmap(mat_scaled, name = "Expression", col = cols, 
                top_annotation = c(ha1.5, ha1, ha2, ha3),
                bottom_annotation = ha,
                show_row_names = FALSE, show_column_names = F)
  
  # plot the heatmap 
  #ppi=300
  pdf(file = paste0("./outputData/WGCNA/ME_vs_trait/ME_heatmaps/heatmap_", module_name, ".pdf"), width=10.5, height=6, pointsize = 1)
  draw(ht1, row_title = "Genes", column_title = paste0(module_name), background = "transparent")
  dev.off()
  
  # Return heatmap
  #draw(ht1, row_title = "Genes", column_title = paste0(module_name))
}
