print("Depicting spatial clustering patterns")
library(Giotto)
library(smfishHmrf)
library("optparse")
python_path = NULL
option_list = list(
  make_option(c("-r", "--resultfile"), type="character", default="results", 
              help="file name for results", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
results_folder=opt$r
if(is.null(python_path)) {
  installGiottoEnvironment()
}
instrs = createGiottoInstructions(show_plot =TRUE,save_plot = TRUE, 
                                  save_dir = results_folder,python_path = python_path)
obj=readRDS(file=paste0(results_folder,"/Giotto_obj.rds"))
obj <- createSpatialNetwork(gobject = obj, method = 'kNN', k = 5, name = 'spatial_network')

my_spatial_genes = readRDS(file=paste0(results_folder,"/spatial_genes.rds"))

#HMRF
hmrf_folder = fs::path(results_folder,'11_HMRF')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
# do HMRF with different betas
HMRF_spatial_genes = doHMRF(gobject = obj, expression_values = 'scaled',spatial_genes = my_spatial_genes,spatial_network_name = 'Delaunay_network',k = 5,betas = c(28), output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k5_scaled'))

obj = addHMRF(gobject = obj,HMRFoutput = HMRF_spatial_genes,k = 5, betas_to_add = c(28),hmrf_name = 'HMRF_2')
spatPlot(gobject = obj, cell_color = 'HMRF_2_k5_b.28', point_size = 3, coord_fix_ratio = 1, save_param = c(save_name = '11_HMRF_2_k5_b.28', base_height = 3, base_width = 9, save_format = 'png'),save_plot=TRUs)
saveRDS(obj,file=paste0(results_folder,"/Giotto_obj.rds"))
