getDF = function(gene1, gene2 = NULL) {
  if (length(gene1) > 1) {
    gene2 = gene1[2]
    gene1 = gene1[1]
  }
  
  gpair = paste0(gene1,"_",gene2)
  if (!gpair %in% rownames(wcors_all)) {
    gpair = paste0(gene2,"_",gene1)
  }
  
  wcor = wcors_all[gpair,]
  
  df_res = data.frame(x = coords[,"x"], 
                      y = coords[,"y"], 
                      g1 = expr[gene1,], 
                      g2 = expr[gene2,],
                      wcor = wcor, 
                      W_min = W[which.min(wcor),], 
                      W_max = W[which.max(wcor),])
  return(df_res)
}
basicplotFunction = function(gene1, gene2 = NULL) {
  
  require(ggforce)
  require(patchwork)
  require(ggpubr)
  
  if (length(gene1) > 1) {
    gene2 = gene1[2]
    gene1 = gene1[1]
  }
  
  df_res = getDF(gene1, gene2)
  
  t = theme(legend.key.width = unit(0.5, "inches")) +
    theme(plot.title = element_text(size = 20)) +
    theme(axis.title = element_text(size = 15))
  
  g_gene1 = ggplot(df_res, aes(x = -x, y = y)) + 
    geom_point(aes(colour = g1), size = 5) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle(gene1) +
    labs(colour = "") +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    
    scale_color_viridis_c(breaks = c(0,max(df_res$g1)),
                          limits = c(0,max(df_res$g1)),
                          labels = c("Low","High")) +
    
    t +
    coord_fixed() +
    guides(colour = guide_colourbar(title.position = "top",
                                    title.hjust = 0.5)) +
    theme(legend.title=element_text(size=15)) +
    labs(colour = "Expression") +
    NULL
  
  g_gene2 = ggplot(df_res, aes(x = -x, y = y)) + 
    geom_point(aes(colour = g2), size = 5) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle(gene2) +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    labs(colour = "") +
    theme(legend.position = "bottom") +
    scale_color_viridis_c(breaks = c(0,max(df_res$g1)),
                          limits = c(0,max(df_res$g1)),
                          labels = c("Low","High")) +
    t +
    coord_fixed() +
    theme(legend.position = "none") +
    NULL
  
  g2 = ggplot(df_res, aes(x = -x, y = y, fill = wcor)) + 
    geom_voronoi_tile(max.radius = 1) +
    theme_minimal() + 
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    labs(colour = "",fill = "") +
    labs(colour = "Local correlation",fill = "Local correlation") +
    ylab("") +
    xlab("") +
    ggtitle("Correlation of both genes") +
    scale_alpha_continuous(range = c(0,0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) + 
    t +
    coord_fixed() +
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = 0.5)) +
    theme(legend.title=element_text(size=15)) +
    NULL
  g_gene1_leg = as_ggplot(get_legend(g_gene1))
  g2_leg = as_ggplot(get_legend(g2))
  
  scater::multiplot(g_gene1 + theme(legend.position = "none") + 
                      theme(plot.margin = margin(10,0,-10,0)),
                    g_gene2 + 
                      theme(plot.margin = margin(10,0,-10,0)),
                    g2 + theme(legend.position = "none") + 
                      theme(plot.margin = margin(10,0,-10,0)),
                    g_gene1_leg,
                    g2_leg,
                    layout = matrix(
                      c(1,1,2,2,3,3,
                        1,1,2,2,3,3,
                        1,1,2,2,3,3,
                        6,4,4,6,5,5), ncol = 6, byrow = TRUE))
}
plotFunction = function(gene1, gene2 = NULL) {
  
  require(ggforce)
  require(patchwork)
  
  if (length(gene1) > 1) {
    gene2 = gene1[2]
    gene1 = gene1[1]
  }
  
  df_res = getDF(gene1, gene2)
  
  g_W_min = ggplot(df_res, aes(x = g1, y = g2)) + 
    geom_point(aes(alpha = W_min, size = W_min), colour = "purple") +
    theme_minimal() +
    xlab(gene1) +
    ylab(gene2) +
    ggtitle("Min") +
    NULL
  g_W_max = ggplot(df_res, aes(x = g1, y = g2)) + 
    geom_point(aes(alpha = W_max, size = W_max), colour = "orange") +
    theme_minimal() +
    xlab(gene1) +
    ylab(gene2) +
    ggtitle("Max") +
    NULL
  
  g_xy = ggplot(df_res, aes(x = x, y = y)) + 
    geom_point(size = 0.1) +
    geom_density_2d(data = subset(df_res, W_max != 0), colour = "orange") +
    geom_density_2d(data = subset(df_res, W_min != 0), colour = "purple") +
    theme_minimal() +
    xlab("x coordinate") +
    ylab("y coordinate") +
    ggtitle("Positions") +
    NULL
  
  g_gene1 = ggplot(df_res, aes(x = x, y = y)) + 
    geom_point(aes(colour = g1), size = 5) +
    
    theme_minimal() +
    ggtitle(gene1) +
    scale_color_gradient2(low = "black", mid = "yellow", high = "red", midpoint = 2) +
    NULL
  
  g_gene2 = ggplot(df_res, aes(x = x, y = y)) + 
    geom_point(aes(colour = g2), size = 5) +
    
    theme_minimal() +
    ggtitle(gene2) +
    scale_color_gradient2(low = "black", mid = "yellow", high = "red", midpoint = 2) +
    NULL
  
  g2 = ggplot(df_res, aes(x = x, y = y, fill = wcor)) + 
    geom_voronoi_tile(max.radius = 1) +
    geom_point(size = 1, colour = "black") +
    theme_minimal() + 
    scale_alpha_continuous(range = c(0,0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1)) + 
    NULL
  
  
  return(g_gene1 + g_gene2 + g2 + 
           g_W_min + g_W_max + g_xy + plot_layout(ncol = 3, nrow = 2, byrow = TRUE))
}
