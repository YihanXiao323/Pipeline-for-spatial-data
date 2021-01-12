library(Giotto)
#load a mel ST data after preprocessing including normalization, dimension reduction, etc.
mel=readRDS(file="Exampledata/melobeject.rds")
#Then test the neighborKMeans
meltest=data.frame(mel@cell_metadata$cell_types)
colnames(meltest)="celltype"
meltest$HMRF=mel@cell_metadata$HMRF_2_k5_b.28
meltest$x=mel@spatial_locs$sdimx
meltest$y=mel@spatial_locs$sdimy
rownames(meltest)=mel@spatial_locs$cell_ID
# calaulate the physical distance for each cell pair
all_dist = dist(meltest[,c(2,3)],method = "maximum")
all_dist = as.matrix(all_dist)

# for each cell, find the nearest neighbours within a specified range#
spot_within = apply(all_dist, 1, function(x){
  #print(quantile(x,0.98))
  return(colnames(all_dist)[order(x)[2:11]])#
})

names(spot_within) = rownames(all_dist)
#cell type aggegation
celltypes=unique(meltest$celltype)
n_celltypes=length(unique(meltest$celltype))
propmatrix=sapply(seq(1,dim(spot_within)[2]),function(x){
  meltest$celltype[spot_within[,x]]
  a=table(meltest[spot_within[,x],"celltype"])
  init=rep(0,n_celltypes)
  names(init)=celltypes
  init[names(a)]=a
  return(init)}
)
propmatrix=t(propmatrix)


# kmeans clustering according to the cell-type proportion of neighbours
kmeans_neighbor_fraction = kmeans(propmatrix, 8, iter.max = 300, nstart = 10)
meltest$kmeanscluster=kmeans_neighbor_fraction$cluster
mel@cell_metadata$kmeans=kmeans_neighbor_fraction$cluster
#plot the KMEANS
spatPlot(gobject = mel, cell_color = 'kmeans', point_size = 3, coord_fix_ratio = 1, save_param = c(save_name = 'kmeans_8', base_height = 3, base_width = 9, save_format = 'png'))
