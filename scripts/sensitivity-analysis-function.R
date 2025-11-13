
sensitivity_analysis_degree <- function(module, no_genes_to_keep){
  #no_genes_to_keep<- 100
  #module <- "2"
  
  message(paste0("running module ", module))
  # get genes of interest to plot 
  genes_to_keep <- select_network_genes(no_genes_to_keep, module)
  # calculate adjacency for genes of interest 
  tom_matrix <- calc_adjacency(genes_to_keep, softThres = 20)
  # loop over several TOM thresholds and visualize degree for each 
  tom_thresh_c <- c(rev(seq(0.5, 0.95, 0.05)))
  
  networks_list <- list()
  degree_tmp <- list()
  for(i in 1:length(tom_thresh_c)){
    networks_list[[i]] <- build_network(module, "#E98113", no_genes_to_keep, tom_thresh_c[[i]])
    degree_tmp[[i]] <- degree(networks_list[[i]])[order(names(degree(networks_list[[i]])))]
  }
  # add names to list for each tom threshold applied 
  names(degree_tmp) <- as.character(tom_thresh_c)
  
  # unpack into data frame 
  degree_df <- do.call(rbind.data.frame, degree_tmp)
  rownames(degree_df) <- names(degree_tmp) 
  colnames(degree_df) <- names(degree_tmp[[2]]) 
  
  # transpose 
  degree_df_t <- as.data.frame(t(degree_df))
  
  # determine top X nodes 
  top_nodes <- names(degree_tmp[[1]][order(degree_tmp[[1]], decreasing = TRUE)])[1:5]
  
  # add genes as own column 
  degree_df_t$genes <- rownames(degree_df_t)
  
  # set colors for plotting 
  degree_df_t$colors <- "gray75"
  degree_df_t$colors[match(top_nodes, rownames(degree_df_t))] <- "red"
  
  degree_df_t$colors[which(degree_df_t$genes==top_nodes[1])] <- "red"
  degree_df_t$colors[which(degree_df_t$genes==top_nodes[2])] <- "cornflowerblue"
  degree_df_t$colors[which(degree_df_t$genes==top_nodes[3])] <- "purple"
  degree_df_t$colors[which(degree_df_t$genes==top_nodes[4])] <- "orange"
  degree_df_t$colors[which(degree_df_t$genes==top_nodes[5])] <- "forestgreen"
  
  
  # set shapes for plotting 
  degree_df_t$shape <- 16
  degree_df_t$shape[which(degree_df_t$genes==top_nodes[1])] <- 17
  degree_df_t$shape[which(degree_df_t$genes==top_nodes[2])] <- 17
  degree_df_t$shape[which(degree_df_t$genes==top_nodes[3])] <- 17
  degree_df_t$shape[which(degree_df_t$genes==top_nodes[4])] <- 17
  degree_df_t$shape[which(degree_df_t$genes==top_nodes[5])] <- 17
  
  degree_df_t$size <- ifelse(degree_df_t$shape==16, 1, 2) 
  
  # convert to long form df 
  data_long <- gather(degree_df_t, tom_thresh, degree, colnames(degree_df_t)[1:10], factor_key=TRUE)
  
  # plot degree vs tom threshold 
  #ppi=300
  pdf(paste0("outputData/WGCNA/network/network_plots/network-degree-sensitivity-analysis-", module, ".pdf"), width=6, height=6, pointsize = 1)
  plot(data_long$degree ~ jitter(as.numeric(as.character(data_long$tom_thresh)), 0.9), 
       pch = data_long$shape, 
       cex = data_long$size,
       xlab = "TOM threshold for edges (percentile)",
       ylab = "Degree",
       main = paste0("Module ", module),
       las=1,
       cex.lab = 1.3, cex.axis = 1.6, cex.main = 1.7,
       col = data_long$colors)
  legend("topright", top_nodes, 
         pch = c(17), 
         col = c("red", "cornflowerblue", "purple", "orange", "forestgreen"), cex=1.05)
  dev.off()
  # return top genes 
  top_nodes
}
