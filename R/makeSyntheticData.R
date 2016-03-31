# Make synthetic lineage-derived single-cell dataset
# Geoffrey F. Schau
# Designs a simulated dataset representing single cells sampled from a heterogeneous tissue 
# with branching lineages

#' makeSyntheticData
#'
#' This function generates a synthetic dataset for use in example documentation
#' @keywords makeSyntheticData
#' @examples
#' makeSyntheticData()

makeSyntheticData <- function(){
  set.seed(1000)
  numCells = 700
  numGroups = 7
  numGenes = numGroups + 1 # +1 for "house-keeping"
  groupNum = floor(numCells/numGroups)
  expressionMatrix = matrix(nrow=numCells, ncol=numGenes, data=0)
  
  # --- Table Labels ---
  colnames(expressionMatrix) = paste('Gene',seq(1,numGenes),sep='')
  rownames(expressionMatrix) = paste('Cell',seq(1,numCells),sep='')
  
  # --- Lineage Diagram ---
  lineages = list(list(2,1),list(3,1),list(4,2,1),list(5,2,1),list(6,3,1),list(7,3,1))
  
  # --- Group-Specific Expression --- 
  for (i in 1:numGroups){ # for each group
    for (j in 1:groupNum){ # for each cell per group
      geneExpression = seq(from=0,to=1,length.out=groupNum)
      expressionMatrix[((i*groupNum)-(groupNum-1)):(i*groupNum),i] = geneExpression
    }
  }
  
  # --- Precusor Expression ---
  for (i in 1:length(lineages)){ # for each lineage
    lineage = unlist(lineages[[i]])
    for (k in 2:length(lineage)){ # for each preceeding part of the lineage
      precLin = lineage[k]
      # NOT RIGHT
      expressionMatrix[((lineage[1]*groupNum)-(groupNum-1)):(lineage[1]*groupNum),precLin] = 1
      #expressionMatrix[((precLin*groupNum)-(groupNum-1)):(precLin*groupNum),lineage[1]] = 1
    }
  }
  
  # --- Add Random Genes ---
  expressionMatrix[,numGenes] = runif(numCells,min=0,max=1)
  
  # --- Save as Table ---
  write.table(expressionMatrix,file='SyntheticData.rda',sep='\t',row.names = TRUE)

}
