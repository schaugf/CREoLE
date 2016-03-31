# Vector-based Single Cell Ordering
# Author: Geoffrey F. Schau and Andrew Adey
# Version 1.3

# Revision History
# v1.1 - July 6, 2015 - Initial Revision
# v1.2 - July 27, 2015 - Code conversion from script to function
# v1.3 - July 29, 2015 - Revised input/output arguments

# Input arguments:
# DM - nx3 matrix containing 3D diffusion map or other low dimension space coordinates for n cells
# k - array of integers specifying which cluster each cell belongs to
# c - array of integers specifying which clusters to order by VbSCO

# Output arguments:
# returnData (data.frame) with the following columns
# - vectorAlign: x,y,and z coordinates of cell projection (columns 1:3)
# - cellNum: cell number from original DM
# - development: each cell's position along the development vector

#' VbSCO
#'
#' The Vector-based Single Cell Ordering (VbSCO) algorithm is designed to project cells of lineage-defining clusters onto centroid-binding vectors
#' @param DM Diffusion Map coordinates
#' @param k K-means cluster array
#' @param c list of cluster that define the lineage
#' @keywords VbSCO
#' @examples
#' VbSCO()

VbSCO <- function(DM,k,c){
  
# Append cluster number to DM
tD = matrix(ncol=3,nrow=nrow(DM),data=0)
tD[1:nrow(DM),1:(min(3,ncol(DM)))] = DM[,1:min(3,ncol(DM))]

DM = data.frame(tD,k,n=1:nrow(tD)) # n indexes sample.ints

# Calculate relevant clusters' centroids
centroids = matrix(ncol=3,nrow=length(c),data=0)

for (i in 1:length(c)){
  cluster = DM[k==c[i],-4]
  centroids[i,1] = sum(cluster[,1])/dim(cluster)[1]
  centroids[i,2] = sum(cluster[,2])/dim(cluster)[1]
  centroids[i,3] = sum(cluster[,3])/dim(cluster)[1]
}

# Select only relevant cell clusters
subDM = DM[k %in% c,]

# ==============================================
# Create planar boundaries 
# ==============================================
# For every cluster that is not the origin or terminal cluster, a decision boundary is needed.
# Decision planes are stored as the normal vector to their respective centroids
# See: http://math.kennesaw.edu/~plaval/math2203/linesplanes.pdf

if (length(c)>=3){ # Two or more vectors connect three or more centroids
  nPlane <- matrix(ncol = 3,nrow=length(c)-2)
  for (i in 1:(length(c)-2)){
    pPlane = centroids[(i+1),] # Passes through centroids, starting at second (first is origin)
    
    # Calculate two unit vectors from each centroid
    v1 = centroids[(i+1),] - centroids[(i),] # center centroid minus preceeding
    v2 = centroids[(i+1),] - centroids[(i+2),] # center centroid minus following
    u1 = v1/sqrt(v1[1]^2 + v1[2]^2 + v1[3]^2) # unit vectors (divide by length)
    u2 = v2/sqrt(v2[1]^2 + v2[2]^2 + v2[3]^2)
    
    # Calculate bisecting vector
    v3 = (u1 + u2) / 2
    
    # Calculate right vector for 3-point plane
    v4 = crossprod(v3,v2)
    
    # Calculate normal vector of plane, which bisects the other two vectors
    nPlane[i,] = crossprod(v4,v3) # will be used to evaluate decision boundary later on
  }
}

# ==============================================
# Vector Allocation and Projection
# ==============================================
# Allocate projection matrix and vector association. v is whichever vector the cell is projected to
# See: http://math.etsu.edu/multicalc/prealpha/Chap1/Chap1-1/printversion.pdf

projMat = matrix(ncol = 3,nrow = nrow(subDM)) # numCells in subset
v = rep(0,nrow(subDM))
subDM = data.frame(subDM,v)

# Allocate cells to projection vector
for (i in 1:nrow(subDM)){
  P = c(subDM[i,1],subDM[i,2],subDM[i,3])
  Pc = subDM$k[i] # cluster
  
  # Check whether cluster is either origin or terminus and allocate accordingly
  if (!(Pc == c[1] | Pc == c[length(c)])){
    # Evaluate decision boundary (if applicable)
    dPlane = nPlane[(which(c==Pc))-1,]
    # Plug in cell location into plane equation and solve (use plane equations from above)
    d = dPlane %*% (P - centroids[(which(c==Pc)),])
    #cat(d)
    if (d>=0){
      subDM$v[i] = which(c==Pc) - 1 # TO SWAP VECTOR ASSIGNMENT
    } else {
      subDM$v[i] = which(c==Pc)
    }
  } else if (Pc == c[1]){
    subDM$v[i] = 1
  } else {
    subDM$v[i] = length(c)-1
  }
  
  # Project each cell onto prescribed vector
  # For a line defined by two points A and B, project point P by:
  # A + dot(AP,AB) / dot(AB,AB) * AB
  A = centroids[(subDM$v[i]),]
  B = centroids[(subDM$v[i]+1),]
  
  # Calculate necessary projection vectors
  AP = c(A[1]-P[1],A[2]-P[2],A[3]-P[3])
  AB = c(A[1]-B[1],A[2]-B[2],A[3]-B[3])
  
  # Calculate projection and append clusters
  projMat[i,] = A + ((AP%*%AB)/(AB%*%AB)) * -AB
}

# I don't think this is right...
projMat = data.frame(projMat,k=subDM$k,n=subDM$n,v=subDM$v)

# ==============================================
# Organize Development Vector
# ==============================================
# Begin origin at centroid. Remove cells "behind" origin by planar boundary
Ori = centroids[1,] - centroids[2,]
dd = numeric(nrow(projMat))
for (i in 1:nrow(projMat)){
  if (projMat[i,4] == c[1]){
    dd[i] = Ori %*% (as.double(projMat[i,1:3]) - centroids[1,])
  }
}
# Filter out cells behind origin centroid
projMat = projMat[dd<=0,]

# Calculate distances between each of the centroids following the origin
# Cells on vector 2 need add length of vector 1; cells on vector 3 need to add length of vector 1 and 2...
centDiff = diff(centroids)
clen = numeric(dim(centDiff)[1])
if (length(clen)!=1){
  for (i in 1:(length(clen)-1)){
    # We start at two because we want the first element to be zero
    # And we don't consider the last element in centDiff because the vector
    # goes to infinity without making another "turn"
    clen[i+1] = sqrt(centDiff[(i),1]^2 + centDiff[(i),2]^2 + centDiff[(i),3]^2)
  }
}

# Measure each cell's position along the development vector, d
d = numeric(nrow(projMat))
for (i in 1:nrow(projMat)){
  # Determine point of origin for each cell's vector (previous centroid)
  # For any given vector, add the total length of the previous vector(s) to it
  # Of course, add zero for first vector (clen[1])
  # For a given vector, measure distance relative from preceding centroid
  # So vector 1 measures from centroid 1, vector 2 measures from centroid 2
  # In this case, the range of subDM$v should be 1 less than dim(centroids)[1]
  cellCentroid = centroids[projMat$v[i],]
  # Find the vector from the cell to that centroid
  cellCentroidVector = projMat[i,1:3] - cellCentroid
  distFromCentroid = sqrt(cellCentroidVector[1]^2 + cellCentroidVector[2]^2 + cellCentroidVector[3]^2)
  # Add the lengths of preceding vectors to distance from centroid
  d[i] = distFromCentroid + sum(clen[1:projMat$v[i]])
}

# Normalize development vector from 0 to 1
d = (d-min(d))
d = d / max(d)

# return data frame and order cells by position along development vector
returnData = data.frame(vectorAlign = projMat[1:3],cellNum = projMat$n,development = d)
returnData = returnData[order(returnData$development),]

}







