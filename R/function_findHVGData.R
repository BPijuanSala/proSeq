########################################################################################
## Title: methodRNAseq_findHVG.R
## Author: Blanca Pijuan-Sala
## Description: Method to find highly variable genes. From Brennecke et al., 2013.
## Date: 30 June 2018
## **proSeq package**
########################################################################################

##----------------------------------------
## findHVG
##----------------------------------------

#' @title Find highly variable genes
#' @description It finds highly variable genes using the normalised data stored in
#' the \code{countsNorm} slot within the \code{RNAseq} object provided. It follows
#' Brennecke et al., Nature Methods, 2013. object@SpikeIn vector is required. 
#' If any entry is true, it will take spike-ins to calculate the fit. If none is true,
#' it will take biological genes to calculate the fit.
#' @param data \code{RNAseq} object.
#' @param signThres Threshold to adjust for multiple testing with the Benjamini-Hochberg 
#' method. Default: 0.1 (cut at 10 percent).
#' @param outputPlots state directory where you want to output the 
#' plots. Default: Current directory.
#' @param plotting It states whether plots are generated or not and you should specify
#' the type. Options: "pdf", "tiff", "no". Default: "pdf". Default: TRUE
#' @param colVarGenes Color for variable genes. Default: deeppink.
#' @param minQuantCv2 Lower quantile of cv2 values to discard to fit the data. Default: 0.2. This means that those
#' genes with cv2 (across cells) in the lower 0.2 quantile will be excluded for the fit.
 #' @param maxQuantCv2 Upper quantile of cv2 values to discard to fit the data. Default: 0.8. This means that those
#' genes with cv2 (across cells) in the top 0.8 quantile will be discarded for the fit.
#' @param minQuantMeans Lower quantile of means to discard to fit the data. Default: 0.2. This means that those
#' genes with mean (across cells) in the 0.2 lower quantile will be discarded for the fit.
#' @param maxQuantMeans Upper quantile of means to discard to fit the data. Default: 0.8. This means that those
#' genes with mean (across cells) in the 0.8 upper quantile will be discarded for the fit.
#' @return  vector of highly variable genes. This should be stored in \code{genesHVG}
#' slot within \code{RNAseq} object.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname findHVGData
#' @importFrom matrixStats rowVars
#' @importFrom statmod glmgam.fit
#' @importFrom graphics lines plot par 
#' @importFrom grDevices pdf tiff dev.off
#' @importFrom stats p.adjust pchisq qchisq
#' @export
findHVGData <- function(data, signThres=0.1, outputPlots="./",plotting="pdf",
                   colVarGenes="deeppink",
                   minQuantCv2 = 0.2, maxQuantCv2 = 0.8,
                   minQuantMeans = 0.2,maxQuantMeans = 0.8)  {
            
           
            
            counts <- data
            keep <- rowSums(counts) > 0
            counts <- counts[keep,]
            
            
  
            nSpikes <- counts
            nCounts <- counts
            if (ncol(nCounts)<100){
              sizeFacsSpikes <-  sizeFacsCounts <- scran::computeSumFactors(counts,
                                                                            sizes=seq(5,15,2))
              
            } else {
              sizeFacsSpikes <-  sizeFacsCounts <- scran::computeSumFactors(counts)
                                                                            
            }

            cat("Estimating technical noise with biological genes...\n\n")
            
            
            
            #ESTIMATE TECHNICAL NOISE
            meansSpikes <- rowMeans(nSpikes)
            varsSpikes <- matrixStats::rowVars(nSpikes)
            cv2Spikes <- varsSpikes / meansSpikes^2

            
            #minMeanForFit <- unname( quantile( meansSpikes[ which( cv2Spikes > minQuantCv2 ) ], minQuantMeans ) )
            minMeanForFit <- unname( quantile( meansSpikes, minQuantMeans ) )
            maxMeanForFit <- unname( quantile( meansSpikes, maxQuantMeans ) )
            
            minCVForFit <- unname( quantile( cv2Spikes, minQuantCv2 ) )
            maxCVForFit <- unname( quantile( cv2Spikes, maxQuantCv2 ) )
            
            validmeansSpikes <- unname(which( meansSpikes > minMeanForFit & meansSpikes < maxMeanForFit))
            validCVSpikes <- unname(which( cv2Spikes > minCVForFit & cv2Spikes < maxCVForFit))
            
            useForFit0 <- names(meansSpikes[intersect(validmeansSpikes,validCVSpikes)])
            useForFit <- names(meansSpikes) %in% useForFit0
            
            #useForFit <- meansSpikes >= minMeanForFit
            fit <- statmod::glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansSpikes[useForFit] ), 
                               cv2Spikes[useForFit] ) 
            xi <- mean(1/sizeFacsSpikes)
            a0 <- unname(fit$coefficients["a0"])
            a1 <- unname(fit$coefficients["a1tilde"]-xi)
            
     
            
            cat("Plotting fit...\n\n")
            
            # Prepare the plot (scales, grid, labels, etc.)
            graphics::plot(meansSpikes, cv2Spikes, pch=20, cex=.2, col="blue",log="xy",
                  xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )

            # Plot the fitted curve
            xg <- 10^seq( -2, 6, length.out=1000 )
            graphics::lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=2 )
            # Plot quantile lines around the fit
            df <- ncol(nSpikes) - 1
            graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df, 
                   col="#FF000080", lwd=2, lty="dashed" )
            graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df, 
                   col="#FF000080", lwd=2, lty="dashed" )  
            
            if (plotting == "pdf"){
              grDevices::pdf(paste0(outputPlots, "HVG_technicalNoiseFit.pdf"))
              graphics::par(mar=c(5, 5, 5, 5))
              
              graphics::plot( meansSpikes, cv2Spikes, pch=20, cex=.2, col="blue",
                    log="xy",
                    xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
              # Plot the fitted curve
              xg <- 10^seq( -2, 6, length.out=1000 )
              graphics::lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
              # Plot quantile lines around the fit
              df <- ncol(nSpikes) - 1
              graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df, 
                     col="#FF000080", lwd=2, lty="dashed" )
              graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df, 
                     col="#FF000080", lwd=2, lty="dashed" )  
              dev.off()
            } else if (plotting == "tiff"){
              grDevices::tiff(paste0(outputPlots, "HVG_technicalNoiseFit.tiff"),width = 5.5, height = 5,
                   res = 300, units="in")
              graphics::par(mar=c(5, 5, 5, 5))
              
              graphics::plot( meansSpikes, cv2Spikes, pch=20, cex=.2, col="blue",
                    log="xy", 
                    xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
              # Plot the fitted curve
              xg <- 10^seq( -2, 6, length.out=1000 )
              graphics::lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
              # Plot quantile lines around the fit
              df <- ncol(nSpikes) - 1
              graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df, 
                     col="#FF000080", lwd=2, lty="dashed" )
              graphics::lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df, 
                     col="#FF000080", lwd=2, lty="dashed" ) 
              
              dev.off()
            }
            
            
            cat("Calculating biological variance...\n\n")
            meansGenes <- rowMeans(nCounts)
            varsGenes <- rowVars(nCounts)
            cv2Genes <- varsGenes / meansGenes^2
            
            psia1theta <- mean( 1 / sizeFacsCounts ) + a1 * mean( sizeFacsSpikes / sizeFacsCounts )
            minBiolDisp <- .5^2
            
            m <- ncol(counts)
            cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
            testDenom <- ( meansGenes * psia1theta + meansGenes^2 * cv2th ) / ( 1 + cv2th/m )
            p <- 1 - pchisq( varsGenes * (m-1) / testDenom, m-1 )
            
            padj <- p.adjust( p, "BH" )
            sig <- padj < signThres
            sig[is.na(sig)] <- FALSE
            #table( sig )
            
            nCountsVar <- nCounts[sig, ]
            
            
            cat("Plotting highly variable genes...\n\n")
            
            
            
            plot( meansGenes, cv2Genes,log="xy",
                  pch=20, cex=.3,  col = ifelse( padj < signThres, colVarGenes, "grey" ), 
                  xlab = "Mean normalized read count", ylab = expression(paste("Squared coefficient of variation",(CV^2)) ))
            #axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
            # expression(10^4), expression(10^5) ) )
            #axis( 2, 10^(-2:3), c( "0.01", "0.1", "1", "10" , "100", "1000"), las=2 )
            
            # Plot the plant genes, use a different color if they are highly variable
            #points( meansGenes, cv2Genes, pch=20, cex=.3, 
            #       col = ifelse( padj < .1, "deeppink", "grey" ) )
            
            # Add the technical noise fit, as before
            xg <- 10^seq( -2, 6, length.out=1000 )
            lines( xg, (xi+a1)/xg + a0, col="black", lwd=1 )
            #lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="deeppink", lwd=3 )
          
            
            if (plotting=="pdf"){
              pdf(paste0(outputPlots, "HVG_VariableGenes.pdf"))
              par(mar=c(5, 5, 5, 5))
              
              plot( meansGenes, cv2Genes,log="xy", 
                    pch=20, cex=.3,  col = ifelse( padj < signThres, "deeppink", "grey" ), 
                    xlab = "Mean normalized read count", ylab = expression(paste("Squared coefficient of variation",(CV^2)) ))
              
              # Add the technical noise fit, as before
              xg <- 10^seq( -2, 6, length.out=1000 )
              lines( xg, (xi+a1)/xg + a0, col="black", lwd=1 )
              if (UseSpike==TRUE){
                points( meansSpikes, cv2Spikes, pch=20, cex=.5, col="blue" )
                
              }
                
              
              dev.off()
            } else if (plotting == "tiff"){
              tiff(paste0(outputPlots, "HVG_VariableGenes.tiff"),width = 5.5, height = 5,
                   res = 300, units="in")
              par(mar=c(5, 5, 5, 5))
              
              plot( meansGenes, cv2Genes,log="xy",
                    pch=20, cex=.3,  col = ifelse( padj < signThres, "deeppink", "grey" ), 
                    xlab = "Mean normalized read count", ylab = expression(paste("Squared coefficient of variation",(CV^2)) ))
              
              # Add the technical noise fit, as before
              xg <- 10^seq( -2, 6, length.out=1000 )
              lines( xg, (xi+a1)/xg + a0, col="black", lwd=1 )
              if (UseSpike==TRUE){
                points( meansSpikes, cv2Spikes, pch=20, cex=.5, col="blue" )
                
              }
              
              dev.off()
            }
            

            return(rownames(nCountsVar))
            
           }
         )






