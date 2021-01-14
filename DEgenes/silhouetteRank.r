library(Giotto)
mel=readRDS(file="Exampledata/STmel/melobeject.rds")
spatial_genes<-silhouetteRank(gobject=mel, expression_values="norm",  rbp_p=0.95, examine_top=0.05)
spatGenePlot(mel, expression_values = 'scaled', genes = spatial_genes[1:4]$genes,point_shape = 'border', point_border_stroke = 0.1,show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 2,save_param = list(save_name = 'silRANK_spatialgenes_km'))
