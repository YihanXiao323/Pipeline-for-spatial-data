setwd("~/lab/pipeline")
library(ggplot2)
spotsvalue=apply(obj@raw_exprs,2,sum)
#obj=readRDS(file=paste0(results_folder,"/Giotto_obj.rds"))
obj=readRDS(file="Exampledata/STmel/melobeject.rds")

spotsvalue=apply(obj@raw_exprs,2,sum)
genesvalue=apply(obj@raw_exprs,2,function(x){length(which(x!=0))})

g<-ggplot()
col <- cut(spotsvalue, breaks = c(min(spotsvalue),quantile(spotsvalue,0.25),quantile(spotsvalue,0.5), quantile(spotsvalue,0.75),max(spotsvalue)), 
           labels = c(paste("<",as.character(floor(quantile(spotsvalue,0.25)/1000)),"k",sep = ""),paste(as.character(floor(quantile(spotsvalue,0.25)/1000)),"k-",as.character(floor(quantile(spotsvalue,0.5)/1000)),"k",sep = ""),paste(as.character(floor(quantile(spotsvalue,0.5)/1000)),"k-",as.character(floor(quantile(spotsvalue,0.75)/1000)),"k",sep = ""), 
                      paste(as.character(floor(quantile(spotsvalue,0.75)/1000)),"k-",as.character(floor(max(spotsvalue,0.5)/1000)),"k",sep = "")))
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
