# sigmaOpt.R
# Guassian kernal width (sigma) optimization
# Author: Geoffrey F. Schau
# Version 1.1

# Input Arguments
# data: nxG matrix of cells to be analyzed
# minSig: beginning value (10^sig+stepsize*i)
# maxSig: end value
# stepSize: step size from beginning to end

# Output Arguments
# sigma: optimal sigma value

#' sigmaOpt
#'
#' sigmaOpt is a course optimization routine for gaussian kernal width selection for use in calculating diffusion map coordinates
#' @param data scaled data as an n x G data frame object
#' @param minSigma lower bound for search space
#' @param numSteps number of steps to take in optimizing kernal
#' @param stepSize incremental width of search space
#' @param is.scaled logical value indicates whether data is scaled 
#' @keywords sigmaOpt
#' @examples
#' sigmaOpt()

sigmaOpt <- function(data,minSigma,numSteps,stepSize,is.scaled=FALSE){
  if (!is.scaled){
    scaledData = data.frame(apply(data,2,scale,center=T,scale=T))
  }
  avrdnorm = numeric(numSteps)
  logsigma = numeric(numSteps)
  numCells = dim(data)[1]
  numGenes = dim(data)[2]
  
  distanceMatrix = as.matrix(dist(data,method='euclidean'))
  
  for(i in 1:numSteps){
    #sigma_ = 10^(minSigma+i*(stepSize));
    sigma_ = minSigma+i*(stepSize)
    # Calculate Minkowski distance matrix
    kern = exp(-distanceMatrix/(2*sigma_^2))
    kernSum = colSums(kern)
    # Adapted from MATLAB implementation
    avrdnorm[i] = (sum(log10(kernSum/numCells) / kernSum)) / sum(1/kernSum)
    #logsigma[i] = log10(sigma_)
    logsigma[i] = sigma_
  }
  dim_norm = diff(avrdnorm) / diff(logsigma)
  sigma=10^(logsigma[which.max(dim_norm)])
  return(sigma)
}




