library(Giotto)
mel=readRDS(file="Exampledata/STmel/melobeject.rds")
bin_spatialgenes = binSpect(mel)
spatGenePlot(mel, expression_values = 'scaled', genes = bin_spatialgenes[1:4]$genes,point_shape = 'border', point_border_stroke = 0.1,show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 2,save_param = list(save_name = 'spatialgenes_km'))
