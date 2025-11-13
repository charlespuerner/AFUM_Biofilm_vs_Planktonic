library(ggrepel)
volcano_plot <- function (dat, out_dir, out_name, title, fc_line, alpha, fc_label, alpha_label){
  
  res_df <- as.data.frame(dat)
  res_df$sig <- c()
  res_df$cols <- c()
  for(i in 1:nrow(res_df)){
    if(is.na(res_df$adj.P.Val[i])){
      res_df$sig[i] <- NA
      res_df$cols[i] <- NA
    }
    else if(res_df$adj.P.Val[i]<=alpha & res_df$logFC[i] > fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC >", fc_line)
      res_df$cols[i] <- "indianred"
    } 
    else if(res_df$adj.P.Val[i]<=alpha & res_df$logFC[i] < -fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "indianred"
    } 
    else if(res_df$adj.P.Val[i]<=alpha & res_df$logFC[i]>-fc_line & res_df$logFC[i]<fc_line){
      res_df$sig[i] <- paste0("Q-value <", alpha, ", Log2FC > -", fc_line, "and ", ", Log2FC <", fc_line)
      res_df$cols[i] <- "cornflowerblue"
    } 
    else if(res_df$adj.P.Val[i]>alpha & res_df$logFC[i] > fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC >", fc_line)
      res_df$cols[i] <- "gray47" 
    }
    else if(res_df$adj.P.Val[i]>alpha & res_df$logFC[i] < -fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "gray47" 
    }
    else if(res_df$adj.P.Val[i]>alpha & res_df$logFC[i] < fc_line){
      res_df$sig[i] <- paste0("Q-value >", alpha, ", Log2FC < -", fc_line)
      res_df$cols[i] <- "gray10" 
    }
  }
  
  ppi=300
  pdf(paste0(out_dir, out_name), width=6, height=7, pointsize = 1)
  p = ggplot(res_df, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col=sig),alpha = 0.5, shape = 21, size =2.5, colour = res_df$cols, fill = res_df$cols)  + 
    #xlim(-2.5, 2.5) +
    xlab("Log2 fold change estimate") + ylab("-log10(adj. P-value)") +
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", linewidth = 0.4) + 
    labs(color="") + 
    theme(legend.key = element_blank()) + 
    ggtitle(title) + 
    geom_vline(xintercept = fc_line, colour = "black", linetype="dotted") + 
    geom_vline(xintercept = -fc_line, colour = "black", linetype="dotted") + 
      theme(legend.key.size = unit(1, "cm"), 
          #panel.grid.major = element_line(colour = "#d3d3d3"), 
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6, linetype = "solid"), 
          panel.background = element_blank(),
          axis.text.x=element_text(colour="black", size = 14, hjust = 1), 
          axis.text.y=element_text(colour="black", size = 14), 
          axis.title.x=element_text(colour="black", size = 16), 
          axis.title.y=element_text(colour="black", size = 16), 
          legend.key = element_blank(), 
          legend.text = element_text(colour="black", size = 16), 
          legend.title = element_text(colour="black", size = 16)) + 
    # add labels to genes w/ LFC > 2 and above alpha threshold
    geom_label_repel(data = subset(res_df, logFC > fc_label & adj.P.Val < alpha_label), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = 0.01,
                     nudge_y = 0.01,
                     point.padding = 0.5,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', size = 4) +
    # add labels to genes w/ LFC < -2 and above alpha threshold
    geom_label_repel(data = subset(res_df, logFC < -fc_label & adj.P.Val < alpha_label), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = -0.01,
                     nudge_y = 0.01,
                     point.padding = 0.5,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', size = 4)
    print(p)
    dev.off()
  
}
