make_module_boxplot <- function(module_name, trait, module_eigengenes_df){

  # extract module eigengene and trait info in df 
  module_df <- data.frame("me" = module_eigengenes_df[, module_name], "culture" = meta[, trait])
  
  # generate boxplot 
  boxplot <- ggboxplot(module_df, "culture", "me",
                       color = "culture", palette =c("#F46D43", "#74ADD1"),
                       add = "jitter", shape = "culture", 
                       title = NULL, ylab = "ME Activation Score")+
    theme(
      # plot title
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      # x and y axis titles
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      # axis tick labels
      axis.text.x  = element_text(size = 16, face = "bold"),
      axis.text.y  = element_text(size = 16, face = "bold"),
      # legend title and text (if you have one)
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 14, face = "bold")
    )
  
  # plot the boxplot 
  #ppi=300
  pdf(file = paste0("outputData/WGCNA/ME_vs_trait/ME_vs_trait_boxplots/boxplot_", module_name, ".pdf"), width=4, height=5, pointsize = 1)
  print(boxplot)
  dev.off()
  
  # Return heatmap
  print(boxplot)
}