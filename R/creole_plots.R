# creole_plots.R
# Geoffrey Schau
# This code component is the third of a series of three separate scripts designed to
# collectively encompass the functionality of the CREoLE algorithm.

# Inputs
# 1) creole_map and creole_lines output file

# Outputs
# 1) Boxplots for each gene across each cluster
# 2) Estimated creole trends for each gene in each lineage
# 3) VSCO output for each gene across

# 4) Need to plot the "expected" trends of the synthetic data (plot rawData)
# 4.a) with noise and dropouts

# --- Roxygen ---
#' creole_plots
#'
#' This function calculates consensus trends of gene expression through each lineage established by creole_map
#' @param origin integer
#' @param figW desired figure width
#' @param figH desired figure height
#' @param figU figure measurement units
#' @param figF desired figure format
#' @param figD figure DPI
#' @keywords creole_plot
#' @export
#' @examples
#' creole_plot()
# ---

creole_plots <- function(outDir,figW=86,figH=60,figU='mm',figF ='eps',figD=350){

  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  # --- Check for folder existance --- 
  setwd(outDir)
  load(paste(outDir,'/creole_map_output.RData',sep=''))
  load(paste(outDir,'/creole_lines_output.RData',sep=''))
  
  # check for lineage folders
  figDir = paste(outDir,'/Figures',sep='')
  bD = paste(figDir,'/Boxplots',sep='')
  if (!dir.exists(bD)){
    dir.create(bD)
  }
  for (i in 1:dim(conExp.out)[3]){
    sF = paste(figDir,'/Lineage',i,sep='')
    if (!dir.exists(sF)){
      dir.create(sF)
    }
  }
  
  # --- Boxplots ---
  pData = cbind(data,cluster=as.factor(KM$cluster))
  for (i in 1:length(GoI)){
    cat('Saving Boxplot',i,'\n')
    p = ggplot(pData,aes_string(x='cluster',y=GoI[i])) +
      geom_boxplot(aes(fill=factor(cluster))) +
      ggtitle(GoI[i]) +
      ylab('Expression') +
      xlab('Cluster') +
      theme_minimal() +
      theme(legend.position = "none") 
    print(p) # to console
    # Good, now apply to everything, separate lineage trends, and clean-up VSCO plots. That's it!!
    ggsave(filename=paste(bD,'/',GoI[i],'.',figF,sep=''), plot = last_plot(),width=figW,height=figH,units=figU,dpi=figD)
  }
  
  # --- Lineage Plots ---
  for (i in 1:dim(conExp.out)[3]){ # for each lineage
    for (j in 1:dim(conExp.out)[2]){ # for each gene
      cat('Saving CREoLE Plot Lineage',i,'Gene',j,'\n')
      fN = paste(figDir,'/Lineage',i,sep='')
      p1 = ggplot(as.data.frame(conExp.out[,,i]),aes(x=seq(from=0,to=1,length.out=dim(conExp.out)[1]),y=conExp.out[,j,i])) +
        geom_line() +
        ggtitle(paste('CREoLE: Lineage',i,'Gene',j)) +
        xlab('Development Time') +
        ylab('Expression') +
        expand_limits(y=0) +
        theme_minimal()
      ggsave(filename=paste(fN,'/',colnames(conExp.out)[j],'_creole.',figF,sep=''), plot = p1,width=figW,height=figH,units=figU,dpi=figD)
      
      p2 = ggplot(as.data.frame(VSCO_out[,,i]),aes(x=dev,y=VSCO_out[,(j+1),i])) +
        geom_point() +
        ggtitle(paste('VSCO: Lineage',i,'Gene',j)) +
        xlab('Development Time') +
        ylab('Expression') +
        expand_limits(y=0) +
        theme_minimal()
      ggsave(filename=paste(fN,'/',colnames(conExp.out)[j],'_vsco.',figF,sep=''), plot=p2, width=figW,height=figH,units=figU,dpi=figD)
      #ggsave(filename=paste(fN,'/',colnames(conExp.out)[j],'.pdf',sep=''), plot = grid.arrange(p1,p2,ncol=2))
    }
    
    # --- Composite of all genes in lineage ---
    cat('Composite Plot',i,'\n')
    mDat = cbind(time=seq(from=0,to=1,length.out=dim(conExp.out)[1]),as.data.frame(conExp.out[,,i]))
    mD = melt(mDat,id='time')
    p = ggplot(data=mD, aes(x=time, y=value, group=variable)) + 
      geom_line(aes(color=factor(variable))) +
      ggtitle(paste('Lineage',i)) +
      xlab('Development Time') +
      ylab('Expression') +
      theme_minimal() +
      theme(legend.title=element_blank())
    print(p)
    fN = paste(figDir,'/Lineage',i,sep='')
    ggsave(filename=paste(fN,'/All_Genes.',figF,sep=''), plot=p, width=figW,height=figH,units=figU,dpi=figD)
    #ggsave(filename=paste(fN,'/All_Genes.pdf',sep=''), plot = last_plot())
  }
  
  # --- Gene Plot from each Lineage ---
  gD = paste(figDir,'/Genes',sep='')
  if (!dir.exists(gD)){
    dir.create(gD)
  }
  lLab = paste('Lineage',1:dim(conExp.out)[3])
  for (i in 1:dim(conExp.out)[2]){
    cat('Gene Plot',i,'\n')
    gDat_ = cbind(time=seq(from=0,to=1,length.out=dim(conExp.out)[1]),as.data.frame(conExp.out[,i,]))
    # Each gene independently
    for (j in 2:dim(gDat_)[2]){
      pDF = as.data.frame(cbind(time=gDat_$time,trend=gDat_[,j]))
      g2 = ggplot(data=pDF, aes(x=time, y=trend)) + 
        geom_line() +
        ggtitle(paste('Lineage',j,'Gene',i)) +
        xlab('Development Time') +
        ylab('Expression') +
        scale_color_discrete(labels=lLab) +
        theme_minimal() +
        theme(legend.title=element_blank())
      ggsave(filename=paste(gD,'/',colnames(conExp.out)[i],'lineage',j,'.',figF,sep=''), plot=g2, width=figW,height=figH,units=figU,dpi=figD)
    }
    
    gDat = melt(gDat_,id='time')
    g = ggplot(data=gDat, aes(x=time, y=value, group=variable)) + 
      geom_line(aes(color=factor(variable))) +
      ggtitle(colnames(conExp.out)[i]) +
      xlab('Development Time') +
      ylab('Expression') +
      scale_color_discrete(labels=lLab) +
      theme_minimal() +
      theme(legend.title=element_blank())
    print(g)
    ggsave(filename=paste(gD,'/',colnames(conExp.out)[i],'.',figF,sep=''), plot=g, width=figW,height=figH,units=figU,dpi=figD)
    #ggsave(filename=paste(gD,'/',colnames(conExp.out)[i],'.pdf',sep=''), plot = g)
  }
  
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
      theme_minimal() +
      theme(legend.position="none")
    
    cenA = KM$centers[,i:(i+1)]
    # This isn't right.
    for(j in 1:(nrow(d.mst))){
      #from = c(KM$centers[d.mst$from[j],][i], KM$centers[d.mst$from[j],][(i%%3)+1])
      #to = c(KM$centers[d.mst$to[j],][i], KM$centers[d.mst$to[j],][(i%%3)+1])
      from = c(KM$centers[d.mst$from[j],][i], KM$centers[d.mst$from[j],][i+1])
      to = c(KM$centers[d.mst$to[j],][i], KM$centers[d.mst$to[j],][i+1])
      sD = data.frame(x=from[1],xend=to[1],y=from[2],yend=to[2])
      p[[i]] = p[[i]] + geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=sD)
    }
    ggsave(filename=paste(figDir,'/Projection',i,'.',figF,sep=''), plot = p[[i]], width=figW,height=figH,units=figU,dpi=figD)
    #ggsave(filename=paste(figDir,'/Projection',i,'.pdf',sep=''), plot = p[[i]])
  }
  
  pdf(paste(figDir,'/ClusterLineageMST.pdf',sep=''),width=5,height=5)
  plot(dmst, layout=layout_with_fr, vertex.size=15,vertex.label.dist=0, vertex.color="grey90", edge.arrow.size=0.5)
  dev.off()
  
}









