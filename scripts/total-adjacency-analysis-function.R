# establish function to apply visualization tools to each module 
total_adjacency_degree_analysis <- function(module, softThres, genes_to_label){
  #module<-"2"
  #softThres <- 20
  #genes_to_label <- top_genes_list$ME2
  
  # get genes of interest to plot 
  genes <- names(net$colors)[which(as.character(net$colors) == module)]
  # calculate adjacency for genes of interest 
  adj <- adjacency(norm_counts[, colnames(norm_counts) %in% genes], power = softThres , type = "signed") 
  adj_sums = rowSums(adj)
  adj_df <- data.frame("gene" = names(adj_sums), total_adjacency = adj_sums)
  adj_df <- adj_df[order(adj_df$total_adjacency, decreasing=TRUE),]
  
  adj_df$shape <- 1
  adj_df$shape[match(genes_to_label, adj_df$gene)] <- 17
  adj_df$size <- 1.1
  adj_df$size[match(genes_to_label, adj_df$gene)] <- 2
  
  
  adj_df$colors <- "gray75"
  adj_df$colors[match(genes_to_label, adj_df$gene)[1]] <- "red"
  adj_df$colors[match(genes_to_label, adj_df$gene)[2]] <- "cornflowerblue"
  adj_df$colors[match(genes_to_label, adj_df$gene)[3]] <- "purple"
  adj_df$colors[match(genes_to_label, adj_df$gene)[4]] <- "orange"
  adj_df$colors[match(genes_to_label, adj_df$gene)[5]] <- "forestgreen"
  
  #ppi=300
  pdf(paste0("outputData/WGCNA/network/network_plots/total-adjacency-analysis-", module, ".pdf"), width=5, height=5, pointsize = 1)
  plot(adj_df$total_adjacency, 
       pch = adj_df$shape, 
       cex = adj_df$size,
       xlab = "Ordered nodes",
       ylab = "Total within-module adjacency",
       main = paste0("Total adjacency - ", "Module ", module),
       las=1,
       cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5,
       col = adj_df$colors)
  legend("topright", genes_to_label, 
         pch = c(17), 
         col = c("red", "cornflowerblue", "purple", "orange", "forestgreen"), cex=1.08)
  dev.off()
  
  write.csv(adj_df, paste0("outputData/WGCNA/network/networkSummarySheets/networktotal-adjacency-analysis-", module, ".csv"))
  
  as.data.frame(adj_df)
}
