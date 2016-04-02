# creole_lines
# Geoffrey Schau
# This code component is the second of a series of three separate scripts designed to
# collectively encompass the functionality of the CREoLE algorithm.

# Inputs
# 1) creole_map output file
# 2) MST
# 3) Origin

# Outputs
# 1) NxGxL matrix of N (user defined) consensus resolution, G genes, and L lineages
# 2) List of enriched genes for plotting

# Staging
# 1) Initialization
#   - input parameters
#   - set working directory
#   - load requisite packages and functions
# 2) Define lineages
# 3) Compute consensus trends
#   - Compute non-consensus trends for comparison
# 4) Identify enriched genes
#   - Cis enrichment
#   - Trans enrichment
# 5) Return output to console

# To Do
# - Bootstrapping requires sampling WITH replacement
# - option to define origin or select putative origin
# - option for q-value cutoff
# - q-value/bonferroni multiple testing
# - Returns list of enriched genes
# - Lineage enrichment to identify set of enriched genes
# - Cis/trans branch specifics
# - Return putative origin, option to choose origin in creole_lines
# - Generates an enrichment table of gene, cluster0, cluster1, q-value, and fold change (*-1/n) if n<1
# - Return union set of genes with q<significance

# --- Roxygen ---
#' creole_lines
#'
#' This function calculates consensus trends of gene expression through each lineage established by creole_map
#' @param origin integer value corresponding to the origin cluster of the data. Defaults to putative origin identified in creole_map
#' @param subFrac sub-fraction of cells selected at each iteration of expression trend estimation, given as percent. Defaults to 0.66
#' @param numPts number of points or resolution of final output trend estimation. Defaults to 10000
#' @param numIters number of estimation iterations. Defaults to 100
#' @param GoI List of Genes of Interest. 'NA' results in all genes, which could significantly increase processing time
#' @keywords creole_lines
#' @export
#' @examples
#' creole_lines()
# ---

creole_lines <- function(outDir,origin=NA,subFrac=0.6,numPts=10000,numIters=100,GoI=NA){
# 
#   subFrac=0.7
#   numPts=10000
#   origin=9
#   numIters=10
#   GoI=NA
  
  cat('Initializing...')
  library(ggplot2)
  library(data.table) 
  library(igraph) 
  library(rgl)
  library(ggplot2)
  library(gridExtra)
  library(splines)

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
  
  # --- Load Data ---
  cat('Loading creole_map data...\n')
  load(paste(outDir,'/creole_map_output.RData',sep=''))
  rawData = data
  if (is.scaled == FALSE){
    scaledData = data.frame(apply(rawData,2,scale,center=T,scale=T))
    scaledData[is.na(colMeans(scaledData))] = 0# Check for empty columns
  } else {
    scaledData = rawData
  }
  if (is.na(origin)){
    origin = putative.origin
  }

  # --- Identify lineages in MST ---
  cat('Calculating Lineages...\n')
  lineages = GetLineages(dmst,origin)

  # --- Calculate Consensus Trends ---
  cat('Calculating Consensus Trends...\n')
  
  if (is.na(GoI)){
    GoI = colnames(rawData)
  }
  
  conExp.out = array(dim=c(numPts,length(GoI),length(lineages)),data=0) # 3D array
  colnames(conExp.out) = GoI
  cData = cbind(KM$cluster,rawData)
  
  for (i in 1:numIters){ # for each iteration
    cat('Iteration',i,'\n')
    sample.ints = sort(sample(1:nrow(cData),floor(nrow(cData)*subFrac)))
    sub.cData = cData[sample.ints,]
    sub.sData = scaledData[sample.ints,]
    sub.k = sub.cData[,1]
    sub.cData = sub.cData[,-1]
    if (dimD == 'DM'){
      sub.DM = DiffMap(sub.sData,sigma)
      } else {
      sub.DM = prcomp(sub.cData)
      sub.DM = sub.DM$x
    }
    for (j in 1:length(lineages)){ # for each lineage
      lineage = lineages[[j]]
      VSCO.out = VbSCO(DM=sub.DM,k=sub.k,c=unlist(lineage))
      VSCO.out$cellNum = sample.ints[VSCO.out$cellNum]
      geneOrder = rawData[VSCO.out$cellNum,]
      for (mm in 1:length(GoI)){ # for each gene
        interp = approx(VSCO.out$development,geneOrder[,GoI[mm]],n=numPts)
        conExp.out[,mm,j] = conExp.out[,mm,j] + interp$y
      }
    }
  }
  
  conExp.out = conExp.out/numIters

  #--- VSCO ---
  cat('Calculating singular trends...\n')
  VSCO_out = array(dim=c(nrow(rawData),length(GoI)+1,length(lineages)),data=NA) # 3D array
  colnames(VSCO_out) = c('dev',GoI)
  
  if (dimD == 'DM'){
    sub.DM = DiffMap(scaledData,sigma)
  } else {
    sub.DM = prcomp(rawData)
    sub.DM = sub.DM$x
  }
  
  for (l in 1:length(lineages)){
    cat('Lineage',l,'\n')
    lineage = lineages[[l]]

    VSCO.out = VbSCO(DM=sub.DM,k=KM$cluster,unlist(lineage))
    geneOrder = data[VSCO.out$cellNum,]
    VSCO_exp = cbind(dev=VSCO.out$development,geneOrder)
    for (mm in 1:ncol(VSCO_out)){ # for each gene
      if (mm == 1){
        VSCO_out[1:nrow(VSCO_exp),1,l] = VSCO_exp[,1]
      } else {
        VSCO_out[1:nrow(VSCO_exp),mm,l] = VSCO_exp[,GoI[mm-1]]
      }
    }
  }
  
  # --- Save Results ---
  cat('Saving Results...')
  save(lineages,conExp.out,VSCO_out,GoI,file='creole_lines_output.RData')
  cat('Done!')
}




