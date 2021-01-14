library(Giotto)
python_path = NULL
results_folder="Exampledata"
if(is.null(python_path)) {
  installGiottoEnvironment()
}
instrs = createGiottoInstructions(show_plot =TRUE,save_plot = TRUE, 
                                  save_dir = results_folder,python_path = python_path)
mel=readRDS(file="Exampledata/STmel/melobeject.rds")

##Nearest neighbor network and clustering
mel <- createNearestNetwork(gobject = mel,dimensions_to_use = 1:15, k = 15)
#Leiden clustering
mel <- doLeidenCluster(gobject = mel, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = mel,cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,save_param = list(save_name = '4_a_UMAP_leiden'))
mel = createSpatialNetwork(gobject = mel, minimum_k = 2, maximum_distance_delaunay = 400)
mel <- createSpatialNetwork(gobject = mel, method = 'kNN', k = 5, name = 'spatial_network')

km_spatialgenes = binSpect(mel)

#HMRF
hmrf_folder = fs::path('Exampledata/11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
my_spatial_genes = km_spatialgenes[1:100]$genes
# do HMRF with different betas
HMRF_spatial_genes = doHMRF(gobject = mel, expression_values = 'scaled',spatial_genes = my_spatial_genes,spatial_network_name = 'Delaunay_network',k = 5,betas = c(28,2,3), output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k5_scaled'))
for(i in seq(28, 32, by = 2)) {
  viewHMRFresults2D(gobject = mel,HMRFoutput = HMRF_spatial_genes,k = 5, betas_to_view = i,point_size = 2)
}# adjust beta
mel = addHMRF(gobject = mel,HMRFoutput = HMRF_spatial_genes,k = 5, betas_to_add = c(28),hmrf_name = 'HMRF_2')
spatPlot(gobject = mel, cell_color = 'HMRF_2_k5_b.28', point_size = 3, coord_fix_ratio = 1, save_param = c(save_name = '11_HMRF_2_k5_b.28', base_height = 3, base_width = 9, save_format = 'png'))
