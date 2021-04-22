#environment
library(Giotto)
library("optparse")
option_list = list(
  make_option(c("-c", "--countsfile"), type="character", default=NULL, 
              help="file name for count-matrix", metavar="character"),
  make_option(c("-l", "--locationfile"), type="character", default=NULL, 
              help="file name for spot-location", metavar="character"),
  make_option(c("-r", "--resultfile"), type="character", default="results", 
             help="file name for results", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$c) & is.null(opt$l)){
  print_help(opt_parser)
  stop("Input files not found", call.=FALSE)
}

expr_matrix<-data.table::fread(file=opt$c)
positions <- data.table::fread(file=opt$l)
head(positions)
results_folder=opt$r

colnames(expr_matrix)[1]<-"genes"
expr_mat = as.matrix(expr_matrix[,-1])
rownames(expr_mat) = expr_matrix$genes
python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

instrs = createGiottoInstructions(show_plot = FALSE,save_plot = TRUE, 
                                  save_dir = results_folder,python_path = python_path)
Giotto_obj <- createGiottoObject(raw_exprs = expr_mat, spatial_locs = positions, instructions=instrs)
#Shows how many genes and cells are lost with combinations of thresholds.
filterCombinations(Giotto_obj, expression_thresholds = c(1, 1), gene_det_in_min_cells = c(20, 20, 20), min_det_genes_per_cell = c(20, 32, 100))
Giotto_obj<-filterGiotto(gobject=Giotto_obj, gene_det_in_min_cells=20, min_det_genes_per_cell=20)
non_mito_genes = grep(pattern = 'mt-', Giotto_obj@gene_ID, value = T, invert = T)
non_mito_or_blood_genes = grep(pattern = 'Hb[ab]', non_mito_genes, value = T, invert = T)
Giotto_obj = subsetGiotto(gobject = Giotto_obj, gene_ids = non_mito_or_blood_genes)

Giotto_obj <- normalizeGiotto(gobject = Giotto_obj, scalefactor = 2000, verbose = T)
Giotto_obj<- addStatistics(gobject = Giotto_obj)
Giotto_obj <- calculateHVG(gobject = Giotto_obj, method = 'cov_groups', zscore_threshold = 0.5, nr_expression_groups = 10, save_param=list(save_name="HVGplot", base_height=5, base_width=5))
gene_metadata = fDataDT(Giotto_obj)
featgenes = gene_metadata[hvg == 'yes']$gene_ID

#PCA
Giotto_obj <- adjustGiottoMatrix(gobject = Giotto_obj, expression_values = c('normalized'), batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),  return_gobject = TRUE, update_slot = c('custom'))
Giotto_obj <- runPCA(gobject = Giotto_obj, expression_values = 'custom', genes_to_use = featgenes, scale_unit = F, center=T, method="factominer")
Giotto::plotPCA(gobject=Giotto_obj,save_param = list(save_name = 'PCA'))
Giotto_obj <- runUMAP(Giotto_obj, dimensions_to_use=1:9, n_components=2)
plotUMAP(gobject=Giotto_obj, point_size=1, save_param = list(save_name = 'UMAP'))
Giotto_obj<-createNearestNetwork(gobject=Giotto_obj, dimensions_to_use=1:9, k=20)
Giotto_obj<-doLeidenCluster(gobject=Giotto_obj, resolution=0.6, n_iterations=10, name="leiden")
plotUMAP(gobject=Giotto_obj, cell_color="leiden", point_size=1, save_param = list(save_name = 'UMAP_leiden'))
markers_scarn=findMarkers_one_vs_all(gobject=Giotto_obj, method="scran", expression_values="custom", cluster_column="leiden", min_genes=5)
markergenes_scran = unique(markers_scarn[, head(.SD, 8), by="cluster"][["genes"]])
plotMetaDataHeatmap(Giotto_obj, expression_values="custom", metadata_cols=c("leiden"), selected_genes=markergenes_scran, save_param = c(save_name = 'metaheatmap_scran'))

#visualize
spatPlot(gobject = Giotto_obj, point_size = 2,save_param = list(save_name = 'spatplot'))

#save
saveRDS(Giotto_obj,file=paste0(results_folder,"/Giotto_obj.rds"))

print("Giotto obeject saved")

library(ggplot2)
obj<-Giotto_obj
spotsvalue=apply(obj@raw_exprs,2,sum)
genesvalue=apply(obj@raw_exprs,2,function(x){length(which(x!=0))})

g<-ggplot()
col <- cut(spotsvalue, breaks = c(min(spotsvalue),quantile(spotsvalue,0.25),quantile(spotsvalue,0.5), quantile(spotsvalue,0.75),max(spotsvalue)), 
           labels = c(paste("<",as.character(floor(quantile(spotsvalue,0.25)/1000)),"k",sep = ""),paste(as.character(floor(quantile(spotsvalue,0.25)/1000)),"k-",as.character(floor(quantile(spotsvalue,0.5)/1000)),"k",sep = ""),paste(as.character(floor(quantile(spotsvalue,0.5)/1000)),"k-",as.character(floor(quantile(spotsvalue,0.75)/1000)),"k",sep = ""),paste(as.character(floor(quantile(spotsvalue,0.75)/1000)),"k-",as.character(floor(max(spotsvalue,0.5)/1000)),"k",sep = "")))
g=g+geom_point(aes(x=obj@spatial_locs$sdimx,y=obj@spatial_locs$sdimy,color=col))
g <- g + labs(x = "x", y = "y", title="Spat plot")
g <- g + scale_colour_hue( name = "Counts per spot")
ggsave(paste(results_folder,"/counts distribution per spot.png",sep = ''),plot = last_plot(), device = NULL)
g

g<-ggplot()
g=g+geom_point(aes(x=obj@spatial_locs$sdimx,y=obj@spatial_locs$sdimy,color=genesvalue))
g <- g + labs(x = "x", y = "y", title="Spat plot")
g <- g + scale_color_gradient(name = "genes detected per spot",low = "cyan",high = "red")
ggsave(paste(results_folder,"/Gene-number distribution per spot.png",sep = ''),plot = last_plot(), device = NULL)
g
