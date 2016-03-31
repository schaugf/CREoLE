# Plot 3D diffusion map with joining mst vectors
# Geoffrey Schau

#' plotDM3D
#'
#' plotDM3D returns opens a new window and a 3-dimensional figure viewer to better visualize the diffusion map and lineage mapping 
#' @param DM Diffusion map coordinates
#' @param KM K-means clustering object
#' @param dmst Minimum Spanning Tree object
#' @param fig.dir Figure directory
#' @keywords plotDM3D
#' 

plotDM3D <- function(init.DM,KM,dmst,fig.dir){
  ndim=3
  open3d()
  numClusters = max(KM$cluster)
  colors = rainbow(numClusters)[KM$cluster]
  #colors = rainbow(nrow(init.DM))
  #plot3d(init.DM[,1:ndim],size=.4,col=colors[KM$cluster],type="s",box=F,labels=F,tick=F,xlab=" ",ylab=" ",zlab=" ",nticks=0)
  plot3d(init.DM,size=.4,col=colors,type="s",box=F,labels=F,tick=F,xlab=" ",ylab=" ",zlab=" ",nticks=0)
  plot3d(KM$centers[,1:ndim],add=T,size=.5,type="s",col="black")
  decorate3d(box=F,xlab=" ",ylab=" ",zlab=" ",axes=F)
  grid3d(c("x","y","z"),at=NULL,col="grey")
  view3d(theta=45,phi=10)
  par3d(windowRect = c(50, 50, 1200, 1200))
  axes3d(edges="bbox",xunit=0,xlen=0,yunit=0,ylen=0,zunit=0,zlen=0)
  
  # Plot joining vector
  d.mst = get.data.frame(dmst)
  for(i in 1:(nrow(d.mst))){
    m=matrix(nrow=2,ncol=3)
    m[1,] = KM$centers[d.mst$from[i],1:ndim]
    m[2,] = KM$centers[d.mst$to[i],1:ndim]
    segments3d(m,col="grey30",lwd=3)
  }
}