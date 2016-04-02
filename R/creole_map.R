# creole_map
# Geoffrey Schau
# This code component is the first of a series of three separate scripts designed to
# collectively encompass the functionality of the CREoLE algorithm.

# Inputs
# 1) Number of Clusters

# Outputs
# 1) 2D projections of the diffusion map with labelled clusters (saved)
# 2) Undirected graphical figure of cluster associations (saved)
# 3) Putative origin (console return)
# 4) Output data file

#' creole_map
#'
#' This function maps input data and calculates associated lineages
#' @param data scaled data as an n x G data frame object
#' @param is.scaled logical input to designate input data as scaled or not (z-transform)
#' @param outDir output directory
#' @param numClusters number of clusters to use for lineage calculation
#' @param numDims number of dimensions to use in reducing data dimensionality
#' @param sigma gaussian kernal width for diffusion calculation
#' @param doplot3d logical whether to launch a window of the 3-D subspace. Defaults to FALSE
#' @param dimD method of dimensionality reduction, currently PCA and Diffusion Mapping are implemented
#' @keywords creole_map
#' @export
#' @examples
#' creole_map()

creole_map <- function(data,is.scaled=FALSE,outDir,numClusters,numDims,sigma,doplot3d=FALSE,dimD='PCA'){
  
  cat('Initializing...\n')
  library(data.table)
  library(igraph)
  library(rgl)
  library(ggplot2)
  library(gridExtra)
  library(splines)
  library(fpc)
  library(cluster)
  library(ica)
  
  # --- Check for folder existance
  cat('Checking filesystem...\n')
  if (!dir.exists(outDir)){
    dir.create(outDir)
  }
  figDir = paste(outDir,'/Figures',sep='')
  if (!dir.exists(figDir)){
    dir.create(figDir)
  }
  setwd(outDir)
  
  # --- Calculate Number of Dimensions ---
  if (missing(numDims)){
    cat('Calculating Dimensionality...\n')
    ndim = function(data){
      percAcct = 0.9 # percent variance accounted for
      ndim = sum(cumsum(prcomp(data)$sdev^2/sum(prcomp(data)$sdev^2)) < percAcct)
      return(ndim)
    }
    numDims = ndim(data)
  } 
  cat('Dimensions:',numDims,'\n')
  
  # --- Generate Initial Reduced Dimensional Mapping ---
  if (dimD == 'DM'){
    # --- Scale Data for Diffusion Mapping
    if (!is.scaled){
      scaledData = data.frame(apply(data,2,scale,center=T,scale=T))
      if (any(is.na(colSums(scaledData)))){
        scaledData = scaledData[,-which(is.na(colSums(scaledData)))]
      }
    }
    if (missing(sigma)){
      cat('Calculating Kernal Width...\n')
      sigma = sigmaOpt(data=scaledData,minSigma=1,numSteps=11,stepSize=1)
    }
    cat('Sigma width:',sigma,'\n')
    # --- Calculate Initial Diffusion Map ---
    cat('Calculating Diffusion Coordinates...\n')
    init.DM = DiffMap(scaledData,sigma)
  } else if (dimD == 'PCA'){ # PCA
    cat('Calculating principal coordinates...\n')
    init.DM = prcomp(data)
    init.DM = init.DM$x
  } else if (dimD == 'ICA'){
    cat('Calculating independent components...\n')
    init.DM = icafast(data,nc = numDims)
    init.DM = init.DM$S
  }

  # --- Calculate Number of Clusters ---
  if (missing(numClusters)){
    cat('Calculating Cluster Density...\n')
    asw = numeric(20)
    for (k in 2:20){
      asw[[k]] = pam(init.DM[,1:numDims], k) $ silinfo $ avg.width
    }
    numClusters = which.max(asw)
  }
  cat('Cluster:',numClusters,'\n')
  
  # --- Cluster Low Dimension Space ---
  cat('Clustering...')
  set.seed(2016)
  KM = kmeans(init.DM[,1:numDims],numClusters) # Initial Clustering of DM

  # --- Calculate Minimum Spanning Tree ---
  cat('Calculating Centroid Tree...\n')
  dmst = makeMST(KM$centers)

  # --- Calculate Putative Origin ---
  putative.origin = max(which(betweenness(dmst)==max(betweenness(dmst))))
  cat('Putative Origin Cluster:',putative.origin,'\n')
  
  # --- Save Data ---
  cat('Saving Data...\n')
  save(data,numClusters,dmst,KM,init.DM,putative.origin,numDims,sigma,is.scaled,dimD,file='creole_map_output.RData')

  # --- Plot Results ---
  cat('Plotting Results...\n')
  
  # --- Projections with Vectors ---
  iDM = init.DM[,1:numDims]
  colnames(iDM) = paste('PC',1:numDims,sep='')
  pData = as.data.frame(cbind(iDM, k=KM$cluster))
  p = list()
  d.mst=get.data.frame(dmst)
  for (i in 1:(numDims-1)){
    p[[i]] = ggplot(pData,aes_string(x=paste('PC',i,sep=''),y=paste('PC',(i+1),sep=''))) +
      geom_point(aes(color=factor(k))) +
      ggtitle(paste('Projection',i)) +
      theme_minimal()
      
    cenA = KM$centers[,i:(i+1)]
    for(j in 1:(nrow(d.mst))){
      from = c(KM$centers[d.mst$from[j],][i], KM$centers[d.mst$from[j],][(i+1)])
      to = c(KM$centers[d.mst$to[j],][i], KM$centers[d.mst$to[j],][(i+1)])
      sD = data.frame(x=from[1],xend=to[1],y=from[2],yend=to[2])
      p[[i]] = p[[i]] + geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=sD)
    }
    ggsave(filename=paste(figDir,'/Projection',i,'.pdf',sep=''), plot = p[[i]])
  }
  
  # --- Plot MST ---
  plot(dmst, layout=layout_with_fr, vertex.size=15,vertex.label.dist=0, vertex.color="grey90", edge.arrow.size=0.5)
  dev.off()
  
  pdf(paste(figDir,'/ClusterLineageMST.pdf',sep=''))
  plot(dmst, layout=layout_with_fr, vertex.size=15,vertex.label.dist=0, vertex.color="grey90", edge.arrow.size=0.5)
  dev.off()

  # --- Plot 3D ---
  if (doplot3d){
    creole::plotDM3D(init.DM[,1:numDims],KM,dmst,figDir)
  }

  cat('Done!')
}





