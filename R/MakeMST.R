# Generate minimum spanning tree of clusters
# Geoffrey Schau

#' makeMST
#'
#' makeMST calculates a minimum spanning tree of cluster centroids
#' @param centers n x d matrix or dataframe object of n cluster centers in d dimensions
#' @keywords makeMST
#' 

makeMST <- function(centers){
  from=numeric(sum(1:(nrow(centers)-1)))
  to=numeric(sum(1:(nrow(centers)-1)))
  dist=numeric(sum(1:(nrow(centers)-1)))
  idx=1
  for (i in 1:(nrow(centers)-1)){
    for (j in (i+1):nrow(centers)){
      from[idx] = i
      to[idx] = j
      # FIX AND PUT IN LOOP
      sSum=0
      for (k in 1:ncol(centers)){
        sSum = sSum + (centers[i,k]-centers[j,k])^2
      }
      dist[idx] = sqrt(sSum)
      idx = idx + 1
    }
  }
  distMat = data.frame(from,to,dist)
  nidx = 1
  widx = 1
  verts=numeric(2*ncol(distMat))
  W=numeric(ncol(distMat))
  for (i in 1:nrow(distMat)){
    verts[nidx] = distMat$from[i]
    verts[nidx+1] = distMat$to[i]
    nidx = nidx + 2
    W[i] = distMat$dist[i]
  }
  
  dgraph = graph(verts,directed=F)
  E(dgraph)$weight=W
  dmst = as.undirected(mst(dgraph))
  E(dmst)$color = "grey"
  return(dmst)
}