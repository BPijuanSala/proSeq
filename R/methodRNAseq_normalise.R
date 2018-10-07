########################################################################################
## Title: methodRNAseq_normalise.R
## Author: Blanca Pijuan-Sala
## Description: Method to normalise transcriptional single cell data.
## Date: 25 December 2016
## **proSeq package**
########################################################################################


##----------------------------------------
## normalise
##----------------------------------------

#' @title Normalise single cells from \code{RNAseq} object
#' @description It normalises single cells from countsRaw slot within the 
#' \code{RNAseq} object provided. It follows Brennecke et al., Nature Methods, 2013.
#' @param object \code{RNAseq} object.
#' @param genes vector of genes to normalise. Default: NULL (It will take all the genes of
#' the \code{countsRaw} slot).
#' @param cells vector of cells to normalise. Default: NULL. This will first look at the
#' \code{cellsFailQC} slot and if empty, it will take all the cells in the 
#' \code{countsRaw} slot
#' @param genes.SizF Genes to use to compute genes size factors for each cell. Default: Vector 
#' of "genes" that are not ERCCs.
#' @param ercc.SizF ERCCs to use to compute ERCC size factors for each cell. Default: Vector 
#' of "erccs".
#' @return  List with normalised counts Matrix (countsNorm) and size factors for biological
#' genes (sizeFactorsGenes) and spike-ins (sizeFactorsSpikeIn). countsNorm should be stored in 
#' \code{object@countsNorm} slot and the size factors should be stored in 
#' \code{object@CellsSizeFac} slot as a list (with the two types of size factors).
#' @author Blanca Pijuan Sala.
#' @export
#' @importFrom scran computeSumFactors
#' @importFrom stats cor
#' @rdname normalise
setGeneric("normalise", function(object,cells=NULL,genes=NULL,
                                 genes.SizF=NULL,ercc.SizF = NULL
                                 )
  standardGeneric("normalise"))

#' @rdname normalise
#' @export
setMethod("normalise", "RNAseq",
          function(object, cells=NULL,genes=NULL,
                   genes.SizF=NULL,ercc.SizF = NULL) {

            cat("Setting up variables...\n\n")
            
            if (length(object@SpikeIn) == 0){
              cat("Please provide a logical vector (TRUE/FALSE) and named with cellNames
                  stating whether the cell is spike-in (TRUE) or not (FALSE). Place it in
                  the RNAseq object 'SpikeIn' slot. If no gene is spike-in, set all 
                  cells to FALSE. Spike-in conditions will not be taken into account.
                  \n\n")
              
            }
            if(is.null(cells)){
              if (length(object@cellsFailQC)==0){
                cellNames <- colnames(object@countsRaw)
                
              } else { 
                cellNames <- names(object@cellsFailQC)[which(object@cellsFailQC==FALSE)]
                
                }
            
            } else {
              cellNames <- cells
            }
              
            if (is.null(genes)){
              geneNames <- rownames(object@countsRaw)
              
            } else {
              geneNames <- genes
            }

            # Get number of genes
            rawCounts <- object@countsRaw[geneNames,cellNames]
            
            spikeins <- object@SpikeIn[which(names(object@SpikeIn) %in% geneNames)]
            counts <-rawCounts[which(spikeins == FALSE),]
            erccs <-rawCounts[which(spikeins == TRUE),]
            

            if (is.null(genes.SizF)){
              genes.SizF <- rownames(counts)
            }
            if (is.null(ercc.SizF)){
              ercc.SizF <- rownames(erccs)
            }
            
            

            cat("Computing size factors...\n\n")
            sizeFacs <- scran::computeSumFactors(counts[genes.SizF,])
            names(sizeFacs) <- colnames(counts)
          
            cat("Plotting size factors against library size (should be correlated) ...\n\n")
            #I always plot the deconvolution size factors against the 
            #library sizes, just as a sanity check. They should be correlated
            correl = round(stats::cor(sizeFacs,colSums(data.frame(counts[genes.SizF,]))),2)
            plot(x=sizeFacs,y=colSums(data.frame(counts[genes.SizF,])),ylab="library size",
                                      xlab="size factor",pch=20,main=paste("correlation:",correl))
            nCounts <- t( t(counts) / sizeFacs)
            
            if (any(spikeins) == TRUE){
              
              erccs_sum <- colSums(erccs)
              if (length(ncol(erccs) > 0)){
                sfsERCC <- scran::computeSumFactors(rawCounts[ercc.SizF,])
                
                sizeFacsSpikes <- sfsERCC
                names(sizeFacsSpikes) <- colnames(counts)
                
                #normalised ERCCs
                nERCCs <- t(t(erccs)/sizeFacsSpikes)
                normCounts <- rbind(nCounts,nERCCs)
                
              } else {
                sizeFacsSpikes <- NULL
                normCounts <- nCounts
                
                
              }
 
              
              
              normCounts <- normCounts[geneNames,cellNames]
              
              cat("DONE!\n\n")
              
              return(list(
                countsNorm = normCounts,
                sizeFactorsGenes = sizeFacs,
                sizeFactorsSpikeIn = sizeFacsSpikes 
              ))
            } else {
              normCounts <- nCounts
              
              normCounts <- normCounts[geneNames,cellNames]
              
              cat("DONE!\n\n")
              
              return(list(
                countsNorm = normCounts,
                sizeFactorsGenes = sizeFacs

              ))
            }

            
              
              }
              )






