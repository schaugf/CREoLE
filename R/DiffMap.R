# DiffMap.R
# Diffusion Map Generation from Single-Cell Data
# Author: Geoffrey F. Schau
# Version 1.1

# Input Arguments
# data: nxG matrix of cells to be analyzed

# Output Arguments
# diffMap: diffusion coordinates of eigen vectors (decreasing in size)

#' DiffMap
#'
#' DiffMap returns the diffusion map corrdinates calculated from input data and sigma, the gaussian kernel width
#' @param data scaled data as an n x G data frame object
#' @param sigma kernel width
#' @keywords Diffusion Map

DiffMap <- function(data,sigma){
  
  numCells = dim(data)[1]
  numGenes = dim(data)[2]
  
  distanceMatrix=as.matrix(dist(data),method='minkowski')
  
  if (missing(sigma)){
    sigma = mean(apply(distanceMatrix,2,sd))/ncol(distanceMatrix)  
  }
  
  # Resolve Guassian kernal
  kern = exp(-distanceMatrix/(2*sigma^2))
  # Normalization
  kernSum = colSums(kern)
  kernSq = outer(kernSum,kernSum,"*")^0
  diag(kern) = 0
  normDistance = kern / kernSq
  distanceSum = diag(1/colSums(normDistance)) %*% normDistance
  
  distanceSum[which(is.na(distanceSum))]=0
  distanceSum[which(is.infinite(distanceSum))]=0
  
  # Calculate eigen values and vectors
  eigens = eigen(distanceSum)
  #eigens = eigen(normDistance)
  
  # Sort and return three largest eigen vectors
  eigenVals = matrix(ncol=numCells,nrow=numCells,data=0)
  diag(eigenVals) = eigens$values
  sortedVals = sort.int(Re(diag(eigenVals)),index.return=T,decreasing=T)
  diffReturn = eigens$vectors[,sortedVals$ix]
  diffReturn = diffReturn[,-1]
} 

