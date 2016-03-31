# Idenfity lineages given an mst and origin
# Geoffrey Schau

#' GetLineages
#'
#' This function identifies the lineages from a given MST and cluster of origin. Terminal nodes are defines as all nodes not the origin with a singular edge. getLineages returns all possible lineages beginning at the origin
#' @param dmst Minimum Spanning Tree of cluster centroids objec
#' @param origin.cluster Cluster of origin 
#' @keywords getLineages
#' @examples
#' GetLineages()

GetLineages <- function(dmst,origin.cluster){
  terminal.clusters = which(degree(dmst)==1)
  terminal.clusters = setdiff(terminal.clusters,origin.cluster)
  #all.paths = list(numeric(length(terminal.clusters)))
  all.paths = length(terminal.clusters)
  
  for (i in 1:length(terminal.clusters)){
    path.obj = get.shortest.paths(dmst,from=origin.cluster,to=terminal.clusters[i])
    #r.path = list(numeric(length(path.obj$vpath[[1]])))
    r.path = numeric(length(path.obj$vpath[[1]]))
    for (j in 1:length(path.obj$vpath[[1]])){
      r.path[j] = path.obj$vpath[[1]][j]
    }
    all.paths[i] = list(r.path)
  }
  return(all.paths)
}