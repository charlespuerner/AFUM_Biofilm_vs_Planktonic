library(ggrepel)
ma_plot <- function (dat, out_dir, out_name, title, alpha, fc_label, alpha_label, ylimit_low, ylimit_up, mean_label){
  
  
  
  #dat <- res
  #alpha <- 0.5
  #fc_label <- 2
  #ylimit_low <- 10
  #ylimit_up <p h
  
  res_df <- as.data.frame(dat)
  res_df$logFC_abs <- abs(res_df$logFC)
  res_df$sig <- c()
  res_df$cols <- c()
  for(i in 1:nrow(res_df)){
    if(is.na(res_df$adj.P.Val[i])){
      res_df$sig[i] <- paste0("Q-value >", alpha)
      res_df$cols[i] <- "gray54"
    }
    else if(res_df$adj.P.Val[i]<=alpha){
      res_df$sig[i] <- paste0("Q-value <", alpha)
      res_df$cols[i] <- "indianred"
    }
    else if(res_df$adj.P.Val[i]>alpha){
      res_df$sig[i] <- paste0("Q-value >", alpha)
      res_df$cols[i] <- "gray54"
    }
  }
  
  ppi=300
  pdf(paste0(out_dir, out_name), width=7.5, height=6, pointsize = 1)
  p = ggplot(res_df, aes(AveExpr, logFC)) +
    geom_point(aes(col=sig),alpha = 0.5, shape = 21, size =2.5, colour = res_df$cols, fill = res_df$cols)  +
    #xlim(1, 5000) +
    ylim(ylimit_low, ylimit_up) +
    ylab("Log2 fold change estimate") + xlab("Normalized mean counts") +
    #geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) +
    labs(color="") +
    theme(legend.key = element_blank()) +
    ggtitle(title) +
    #coord_trans(x="log")+
    #scale_x_continuous(trans='log10') +
    geom_hline(yintercept = 0, colour = "black", linetype="dotted") +
    # add labels to genes w/ LFC > 2 and above alpha threshold
    geom_label_repel(data = subset(res_df, logFC_abs > fc_label & adj.P.Val <= alpha & AveExpr > mean_label), aes(label = gene),
                     min.segment.length = 0,
                     #seed = 42,
                     box.padding   = 0.5,
                     #nudge_x = 0.01,
                     #nudge_y = 5,
                     #point.padding = 0.5,
                     label.size = 0.1,
                     #segment.size = 0.3,
                     #direction = "y",
                     #hjust = "left",
                     segment.color = 'grey50', size = 1.9, max.overlaps = Inf) +
    # add vertical fold change lines
    #geom_vline(xintercept = fc_line, colour = "black", linetype="dotted") +
    #geom_vline(xintercept = -fc_line, colour = "black", linetype="dotted") +
    theme(legend.key.size = unit(1, "cm"),
          #panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"),
          panel.background = element_blank(),
          axis.text.x=element_text(colour="black", size = 14, hjust = 1),
          axis.text.y=element_text(colour="black", size = 14),
          axis.title.x=element_text(colour="black", size = 16),
          axis.title.y=element_text(colour="black", size = 16),
          legend.key = element_blank(),
          legend.text = element_text(colour="black", size = 18),
          legend.title = element_text(colour="black", size = 18))
  print(p)
  dev.off()
  
}
