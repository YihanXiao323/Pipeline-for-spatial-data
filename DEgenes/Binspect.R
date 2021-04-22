library(Giotto)
library("optparse")
option_list = list(
  make_option(c("-r", "--resultfile"), type="character", default="results", 
              help="file name for results", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
results_folder=opt$r

obj=readRDS(file=paste0(results_folder,"/Giotto_obj.rds"))
print("Calculating SpatialDEgenes")
obj = createSpatialNetwork(gobject = obj, minimum_k = 2, maximum_distance_delaunay = 400)
bin_spatialgenes = binSpect(obj)
saveRDS(bin_spatialgenes[1:100]$genes,file=paste0(results_folder,"/spatial_genes.rds"))
spatGenePlot(obj, expression_values = 'scaled', genes = bin_spatialgenes[1:4]$genes,point_shape = 'border', point_border_stroke = 0.1,show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 2,save_param = list(save_name = 'spatialgenes_km'))
