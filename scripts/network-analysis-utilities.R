# establish function to select network genes based on top KME & module 
select_network_genes <- function(n, modules){
  #modules<- c("11")
  #n<-100
  
  kme_tmp <- kme[, c(paste0("MM.", modules))]
  kme_sub <- data.frame(gene = rownames(kme), kme = kme_tmp)
  kme_sub <- kme_sub[kme_sub$gene %in% names(net$colors)[net$colors==modules],]
  kme_sub <- kme_sub[order(kme_sub$kme, decreasing = TRUE),]
  
  kme_sub$gene[1:n]
}

calc_adjacency <- function(genes, softThres){
  #genes <- genes_to_keep
  # calculate adjacency matrix for genes to keep 
  adjacency <- adjacency(norm_counts[, colnames(norm_counts) %in% genes], power = softThres , type = "signed") 
  adjacency[adjacency < 0] = 0
  adjacency[adjacency > 1] = 1
  genes_names <- colnames(as.data.frame(adjacency))
  # calculate TOM 
  TOM = TOMsimilarity(adjacency, TOMType="unsigned")
  adj <- TOM
  colnames(adj) <- genes_names
  adj
}

# 
build_network <- function(module, color_nodes, no_genes_to_keep, tom_thresh_perc){
  #color_nodes <- "#E98113"
  #module <- "5"
  #deg_to_label <- 50
  #tom_thresh_perc <- .95
  #no_genes_to_keep <- 100
  
  # get genes of interest to plot 
  genes_to_keep<- select_network_genes(no_genes_to_keep, module)
  # calculate adjacency for genes of interest 
  tom_matrix <- calc_adjacency(genes_to_keep, softThres = 20)
  # set tom threshold 
  tom_thresh <- quantile(tom_matrix, p = tom_thresh_perc)
  
  # replace all values above threshold with 1
  tom_matrix[tom_matrix > tom_thresh] = 1
  tom_matrix[tom_matrix != 1] = 0
  # draw the network 
  network <- graph.adjacency(tom_matrix)
  network <- igraph::simplify(network)  
  vert <- V(network)
  
  # set color of nodes
  V(network)$color <- color_nodes
  #V(network)$color <- net$colors[names(net$colors) %in% genes_to_keep]
  #V(network)$color <- net$colors[net$colors %in% module]
  # set node names and size
  V(network)$names <- colnames(tom_matrix)
  V(network)$size <- 3
  
  # build df of network features
  kme_sub <- kme[match(colnames(tom_matrix), rownames(kme)), ]
  kme_sub_2 <- as.data.frame(kme_sub[, c(colnames(kme_sub) %in% paste0("MM.", module))])
  colnames(kme_sub_2) <- "module"
  rownames(kme_sub_2) <- rownames(kme_sub)
  
  network
}



plot_network <- function(network, module, deg_to_label, no_genes_to_keep){
  # calculate network degree
  deg <- igraph::degree(network, mode = "all")
  
  # Check the name attribute: use 'name' instead of 'names'
  if("name" %in% vertex_attr_names(network)){
    vertex_names <- V(network)$name
  } else if ("names" %in% vertex_attr_names(network)){
    vertex_names <- V(network)$names
  } else {
    stop("The network graph has no vertex attribute 'name' or 'names'.")
  }
  
  names(deg) <- vertex_names
  
  # remove unconnected nodes
  network <- delete_vertices(network, deg == 0)
  
  # plot network
  set.seed(1234)
  pdf(paste0("outputData/WGCNA/network/network_plots/network-", module,"-95Cutoff", ".pdf"), width=10, height=10)
  
  p2 <- plot(network,
             layout = layout_with_fr(network),
             edge.arrow.size = 0.5,
             edge.arrow.mode = 0,
             vertex.label = ifelse(igraph::degree(network, mode = "all") > deg_to_label, V(network)$name, NA),
             vertex.label.cex = 1,
             vertex.label.color = "black",
             main = paste0("Network degree for ", no_genes_to_keep, " most connected genes (kME) - ", module)
  )
  
  dev.off()
  
  return(p2)
}



# establish function to plot TOM for top genes of each module 
module_heatmap <- function(module_to_vis, no_genes_to_keep, softThres){
  #module_to_vis <- "2"
  #no_genes_to_keep <- 100
  #softThres <- 20
  
  # get genes of interest to plot 
  genes_to_keep <- select_network_genes(no_genes_to_keep, module_to_vis)
  # calculate adjacency for genes of interest 
  adj <- calc_adjacency(genes_to_keep, 20)
  # zero 1 values on diagonal 
  adj = adj - diag(diag(adj))
  
  # set up colors for heatmap 
  cols = colorRamp2(c(0, max(adj)*0.7, max(adj)), c("white", "#D38B83", "#DA2511"))
  
  ha = columnAnnotation(x = anno_text(colnames(adj), 
                                      which="column", 
                                      rot = 45, 
                                      gp = gpar(fontsize = 4)))
  
  # generate heatmap object Mi
  ht1 = Heatmap(adj, name = "TOM", col = cols, 
                #top_annotation = c(ha1, ha2),
                bottom_annotation = ha, 
                show_row_names = FALSE, show_column_names = FALSE)
  
  # plot the heatmap 
  #ppi=300
  pdf(paste0("outputData/WGCNA/network/network_plots/TOM-heatmap-module", module_to_vis, "-top_", no_genes_to_keep, "_genes-.pdf"), width=9, height=9, pointsize = 1)
  draw(ht1, row_title = "Genes", 
       column_title = paste0("Module ", module_to_vis, " - TOM (top ", no_genes_to_keep, " genes)"))
  dev.off()
}










# establish function to apply visualization tools to each module 
net_vis <- function(module_to_vis, no_genes_to_keep, softThres){
  #module_to_vis<-"2"
  #no_genes_to_keep <- 100
  
  # get genes of interest to plot 
  genes_to_keep <- select_network_genes(no_genes_to_keep, module_to_vis)
  # calculate adjacency for genes of interest 
  adj <- calc_adjacency(genes_to_keep, softThres = 20)
  # select threshold 
  tom_thresh <- quantile(adj, p = 0.95)
  # plot network 
  net_me1 <- plot_network(adj, "#E98113", module_to_vis, 50, tom_thresh)
}
