########################################################################################
## Title: function_findHVG_DM.R
## Author: Blanca Pijuan-Sala
## Description: Distance to median method to find highly variable genes.
## Date: 30 September 2017
## **proSeq package**
########################################################################################



#' @title Find HVG with the DM approach
#' @description It finds highly variable genes using the distance-to-median 
#' method proposed in Kolodziejczyk et al. Cell Stem Cell, 2015. I adapted the code from 
#' Ximena Ibarra-Soria.
#' @param data Please provide the data without any transformation (e.g. log10).
#' @param qMinThr Quantile threshold to discard lowly expressed genes.E.g. if 
#' set to 0.9, it will take the data with median above the quantile 0.9 
#' and compute the DM. The ones below (bottom 90\%) will be discarded. Default: 0.9.
#' @param qThr Quantile threshold from which to consider a gene highly variable. E.g. If set 
#' to 0.8, it will take the most HVG genes falling on the top 20\%. Default: 0.8
#' @param plot Logical. Default: FALSE.
#' @return vector of HVG.
#' @author Blanca Pijuan-Sala.
#' @references Kolodziejczyk et al. Cell Stem Cell, 2015.
#' @export
#' @rdname findHVG-DM
#' @importFrom zoo rollapply
#' @importFrom stats quantile median 
#' @importFrom graphics plot abline points

findHVG.DM <- function(data,qMinThr=0.9, qThr=0.8,plot=FALSE){
  counts <- data
  #remove genes that are never expressed
  test<-rowSums(counts)
  remove<-which(test==0)
  if(length(remove)>0) { counts<-counts[-remove,]}
    
  df<-data.frame(mean=apply(counts, 1, function(x) mean(x)), 
                 sd=apply(counts, 1, function(x) stats::sd(x)), 
                 cv2=apply(counts, 1, function(x) (stats::sd(x)/mean(x))^2) )
  data <- cbind(log10(df[,"mean"]), log10(df[,"cv2"]))
  colnames(data) <- c("log.avg","log.cv2")
  row.names(data) <- row.names(df)
  
  roll.median.avg<-zoo::rollapply(data[order(data[,"log.avg"]),"log.avg"], 
                              width=50, FUN=stats::median, by=25)
  roll.median.cv2<-zoo::rollapply(data[order(data[,"log.avg"]),"log.cv2"], 
                             width=50, FUN=stats::median, by=25)
    
  ## filter very lowly expressed genes
  #qMinThr<-0.9
  minMean <- max(roll.median.avg[ roll.median.cv2 >= stats::quantile(roll.median.cv2[roll.median.cv2 != max(roll.median.cv2)], qMinThr) ])
  remove <- which(data[,"log.avg"] < minMean)
  data<-data[-remove,]
  
  ## compute DM
  id <- 1:nrow(data)
  dm.mean <- unlist(sapply(id, function(x){
    exp <- data[x,"log.avg"]
    test <- which(roll.median.avg > exp)
    if(length(test)==0){ window = length(roll.median.cv2) }
    else{ window = min(test) }
    return(data[x,"log.cv2"]-roll.median.cv2[window])
  } ))
  names(dm.mean)<-row.names(data)
  
  # select HVGs
  #qThr <- 0.8
  dm.mean <- dm.mean[order(dm.mean, decreasing = T)]
  high.var <- dm.mean[which(dm.mean >= stats::quantile(dm.mean, qThr))]
    
  # plot
  if(plot){
    par(mfrow=c(1,2))
    graphics::plot(data[,"log.avg"], data[,"log.cv2"], xlab = "mean expression", ylab=expression("CV"^2), pch=16)
    graphics::abline(v=(minMean), col="red", lty=2, lwd=2)
    graphics::points(data[high.var,"log.avg"], data[high.var,"log.cv2"], col="red", pch=16)
    
    plot(data[,"log.avg"], dm.mean, xlab="mean expression", ylab="DM", pch=16)
    points(data[high.var,"log.avg"], dm.mean[high.var], col="red", pch=16)
  }
  return(high.var)
}
