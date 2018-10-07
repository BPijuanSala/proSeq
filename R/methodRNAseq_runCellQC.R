########################################################################################
## Title: methodRNAseq_runCellQC.R
## Author: Blanca Pijuan-Sala
## Description: Method to run QC on transcriptionally-profiled cells.
## Date: 24 December 2016
## **proSeq package**
########################################################################################


##----------------------------------------
## runCellQC
##----------------------------------------

#' @title Perform QC on cells from \code{RNAseq} object
#' @description It performs QC on cells from countsRaw slot within the 
#' \code{RNAseq} object provided. It follows Brennecke et al., Nature Methods, 2013.
#' @param object \code{RNAseq} object.
#' @param annTable Associative matrix or dataframe where column 1 reflects geneIDs
#'  and column2 the associated geneName. Default: Mus musculus geneTable (??geneTable).
#' @param geneInput It states whether the genes in the counts table are
#' IDs ("ID") or associated gene names ("name"). Note that IDs will refer to column
#' 1 of annTable and names to column 2 of annTable. In the case of "names": if the ID 
#' is not found, the gene will be discarded. Default: ID (recommended).
#' @param maxMapSpike Maximum percentage of mapped spikeins. Default: 0.5.If
#'  no gene is spike-in, this condition will not be taken into account.
#'  Check the fMapped2ERCCs plot.
#' @param maxMapmit Maximum percentage mapped mitochondrial reads. Default: 0.1.
#' Check the fMapped2Mit plot.
#' @param minMapreads Minimum number of reads mapping to nuclear genes.
#' Default: 200,000. Check the nNuclearGenes plot. Important: This plot is log-scaled.
#' @param GeneDetection How many reads in a gene are needed to consider
#' a gene detected in a cell. Default: 2. Important for nDetectedGenes plot.
#' @param minGenesExpr Minimum number of genes expressed in one cell at least. 
#' Default: 4000. Check nDetectedGenes plot.
#' @param checkCountsFeat Logical (TRUE/FALSE). Specifies whether the QC should take into account
#' the features of the mapping. Default = TRUE (recommended).
#' @param outputPlots state directory where you want to output the plots.
#' Default: Current directory.
#' @param plotting It states whether plots are generated or not and you should specify
#' the type. Options: "pdf", "tiff", "no". Default: "pdf". 
#' Default: TRUE
#' @param passCol Plotting parameter. Color for cells that have passed QC. 
#' Default: "grey70"
#' @param failCol Plotting parameter. Color for cells that have passed QC. 
#' Default: "red".
#' @param thresCol Plotting parameter. Color for threshold line. Default: "blue".
#' @param lty Plotting parameter. Line shape for threshold. For further information
#' see \code{\link[graphics]{par}}. Default: 3. 
#' @param lwd Plotting parameter. Line width for threshold.For further information
#' see \code{\link[graphics]{par}}. Default: 2.
#' @return vector of logical values where TRUE = fail and FALSE = pass
#' stored in \code{runCellQC} slot from \code{RNAseq} object.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname runCellQC
#' @importFrom graphics plot abline 
#' @importFrom grDevices pdf tiff dev.off
setGeneric("runCellQC", function(object,annTable=geneTable, geneInput="ID",
                                 GeneDetection=2,
                                 minMapreads = 200000, maxMapmit = 0.1, 
                                 maxMapSpike = 0.5, minGenesExpr = 4000,
                                 checkCountsFeat = TRUE,outputPlots ="./",
                                 plotting="pdf",passCol="grey70",failCol="red",
                                 thresCol="blue",lty=3,lwd=2) 
  standardGeneric("runCellQC"))

#' @rdname runCellQC
#' @export
setMethod("runCellQC", "RNAseq",
          function(object,annTable=geneTable, geneInput="ID",
                   GeneDetection=2,
                   minMapreads = 200000, maxMapmit = 0.1, 
                   maxMapSpike = 0.5, minGenesExpr = 4000,
                   checkCountsFeat = TRUE,outputPlots ="./",
                   plotting="pdf",passCol="grey70",failCol="red",
                   thresCol="blue",lty=3,lwd=2) {
            cat(paste("Parameters:","\n"))
            cat(paste("  geneInput:",geneInput,"\n"))
            cat(paste("  GeneDetection:",GeneDetection,"\n"))
            cat(paste("  minMapreads:",minMapreads,"\n"))
            cat(paste("  maxMapmit:",maxMapmit,"\n"))
            cat(paste("  maxMapSpike:",maxMapSpike,"\n"))
            cat(paste("  minGenesExpr:",minGenesExpr,"\n"))
            cat(paste("  checkCountsFeat:",checkCountsFeat,"\n"))
            cat(paste("  outputPlots:",outputPlots,"\n\n"))
            
            cat("Setting up variables...\n\n")
            

            if (length(object@SpikeIn) == 0){
              cat("Please provide a logical vector (TRUE/FALSE) and named with cellNames
                  stating whether the cell is spike-in (TRUE) or not (FALSE). Place it in
                  the RNAseq object 'SpikeIn' slot. If no gene is spike-in, set all 
                  cells to FALSE. Spike-in conditions will not be taken into account.
                  \n\n")
              
            }
            if (checkCountsFeat == TRUE){
              if (nrow(object@CountsFeat) == 0){
                cat("Please provide a counts features dataframe or matrix for in the 
                    CountsFeat slot of the object.\n\n
                    ")
              }
              }
            
            
            
            # Get number of genes
            rawCounts <- object@countsRaw[names(object@SpikeIn),]
            #cat(paste("dimensions rawCounts:",dim(rawCounts)[1],dim(rawCounts)[2],"\n"))
            #cat(paste("dimensions annTable:",dim(annTable)[1],dim(annTable)[2],"\n"))
            # cat(paste("dimensions geneTable:",dim(geneTable)[1],dim(geneTable)[2],"\n"))
            
            
            if (geneInput=="name"){
              geneIDs <- annTable[which(annTable[,2] %in% row.names(rawCounts)),1]
              geneNames <- annTable[which(annTable[,1]==geneIDs),2]
              counts <- countsRaw[geneNames,]
              rownames(counts) <- geneIDs
            } else if (geneInput == "ID"){
              counts <- rawCounts
            } else {
              cat("Please enter 'name' or 'ID' in geneInput parameter.\n\n")
            }
            
            #cat(paste("dimensions counts:",dim(counts)[1],dim(counts)[2]))
            
            cellNames <- colnames(counts)
            
            genes <- rownames(counts)
            numGenes <- length(genes)
            
            
            if (geneInput=="name"){
              htseq.genes <- names(object@SpikeIn)[which(object@SpikeIn==FALSE)]
              htseqIDs <- annTable[which(annTable[,2] %in% htseq.genes),1]
              htseq.data <- counts[htseqIDs,]
              
              
              ercc.genes <- names(object@SpikeIn)[which(object@SpikeIn==TRUE)]
              erccIDs <- annTable[which(annTable[,2] %in% ercc.genes),1]
              
              ercc.data <- counts[erccIDs,]
              
              
            } else if (geneInput == "ID"){
              htseqIDs <- names(object@SpikeIn)[which(object@SpikeIn==FALSE)]
              htseq.data <- counts[htseqIDs,]
              
              erccIDs <- names(object@SpikeIn)[which(object@SpikeIn==TRUE)]
              ercc.data <- counts[erccIDs,]
            } else {
              cat("Please enter 'name' or 'ID' in geneInput parameter.\n\n")
            }
            
            
            #genesInd <- match(c(rownames(htseqQC), rownames(ercc.data)), genes)
            #genes <- genes[-match(c(rownames(htseqQC), rownames(ercc.data)), genes)]
            
            
            ## Get ERCC, QC, mitochondrial and nuclear gene IDs separated out
            
            mitochondrialGenesIDsAnn <- annTable[grep("^mt-", annTable[,2],ignore.case = TRUE),1]
            mitochondrialGenesIDs <- rownames(counts)[which(rownames(counts) %in% mitochondrialGenesIDsAnn)]
            
            nuclearGeneIDsAnn <- annTable[-match(mitochondrialGenesIDs, annTable[,1]),1]
            nuclearGeneIDs <- rownames(counts)[which(rownames(counts) %in% nuclearGeneIDsAnn)]
            
            
            if (checkCountsFeat == TRUE){
              nNoFeatureReads <- object@CountsFeat[1,colnames(counts)] # From HT-Seq - number of no feature reads
              nUnalignedReads <- object@CountsFeat[4,colnames(counts)] # From HT-Seq - number of unaligned reads
              nAmbigReads <- object@CountsFeat[2,colnames(counts)] # From HT-Seq -  number of ambiguous reads
              nLowQuality <- object@CountsFeat[3,colnames(counts)] # From the HT-Seq
            }
            
            cat("Computing QC...\n\n")
            
            
            ## QC.
            if (checkCountsFeat == TRUE){
              TotalReads <- rbind(counts,object@CountsFeat[,colnames(counts)])
              nTotalReads <- colSums(TotalReads) # Number of total reads
            } else {
              nTotalReads <- colSums(counts) # Number of total reads
              
            }
            nMappedReads <- colSums(counts[c(nuclearGeneIDs,erccIDs,mitochondrialGenesIDs),]) # Number of mapped reads
            nGenes <- colSums(counts[c(nuclearGeneIDs, mitochondrialGenesIDs),]) # Number of reads mapping to genes/mitchondrial genes
            
            
            nMitochondrialGenes <- colSums(counts[mitochondrialGenesIDs,]) # Number of reads mapping to mitochondrial genes
            nERCC <- colSums(counts[erccIDs,]) # Number of reads mapping to ERCCs
            nNuclearGenes <- colSums(counts[c(nuclearGeneIDs),]) # Number of reads mapping to nuclear genes
            
            #rpm <- t(t(1000000 * counts[nuclearGeneIDs,]) / nNuclearGenes) # Convert into reads per million?
            #nLowCoverageReads <- apply(rpm, 2, function(x) sum(x > rpmThresh)) # Number of genes per cell where rpm > 10
            
            nDetectedGenes <- colSums(counts[c(nuclearGeneIDs), ] >= GeneDetection) # Number of genes with more than GeneDetection reads.
            
            
            # Add thresholds for quality control
            failQC <- names(which(nNuclearGenes < minMapreads)) # Need at least 200,000 reads mapping to nuclear genes
            failQC <- c(failQC, names(which(nMitochondrialGenes/nMappedReads >= maxMapmit))) # Want <10% mapped reads mitochondrial
            if (sum(which(object@SpikeIn ==TRUE))>0){
              failQC <- c(failQC, names(which(nERCC/nMappedReads > maxMapSpike))) # No more than half mapped reads ERCC
            }
            
            failQC <- c(failQC, names(which(nDetectedGenes < minGenesExpr))) # Want at least 4000 genes expressed

            failQC <- unique(failQC)
            
            
            # How many cells fail QC ------------------------------------------
            failQCBool <- cellNames %in% failQC
            names(failQCBool) <- cellNames
            cat(paste("Number cells failed:", sum(failQCBool),"\n\n"))
            
            cat("Generating plots...\n\n")
            # Generate quality control plots
            plotList <- list("fMappedReads", 
                             "fMapped2Nuclear", 
                             "fMapped2ERCC",
                             "fMapped2Mit", 
                             "nTotalReads", 
                             "nNuclearGenes", 
                             "nDetectedGenes")
            plotNames<- c("Mapped reads/Total reads",
                          "Reads mapping to nuclear genes/Mapped reads", 
                          "Reads mapping to Spike-in/Mapped reads", 
                          "Reads mapping to mitochondrial/Mapped reads", 
                          "log10 Total reads", 
                          "log10 Number of reads mapping to nuclear genes", 
                          "Number of detected genes")
            
            logList <- c("", "", "", "", "y", "y", "")
            
            thresholdList <- list(F, F, maxMapSpike, maxMapmit, F, minMapreads, minGenesExpr)
            
            df_qc <- data.frame(nMappedReads/nTotalReads, 
                                nNuclearGenes/nMappedReads, 
                                nERCC/nMappedReads,
                                nMitochondrialGenes/nMappedReads, 
                                nTotalReads, 
                                nNuclearGenes, 
                                nDetectedGenes, 
                                1:length(nDetectedGenes))
            colnames(df_qc) <- c(plotList, "CellIndex")
            
            #plotCols <- rep(c("grey70", "grey20"), each = 96, len = 1920)
            plotCols <- rep(c(passCol), len = ncol(counts))
            names(plotCols)<-cellNames
            plotCols[failQC] <- failCol
            plotCols <- sort(plotCols,decreasing=FALSE)
            
            #pretty_plot <- function(){
            #  theme_bw()+
            #    theme(panel.border = element_rect(colour = "black", size=1), panel.grid.major = element_blank(),
            #          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
            #}
            
            
            for (i in 1:length(plotList)){
              
              if (logList[i] == "y") {
                y_i <- log(df_qc[names(plotCols),plotList[[i]]])
                thres_i <- log(thresholdList[[i]])
              } else {
                y_i <- df_qc[names(plotCols),plotList[[i]]]
              thres_i <- thresholdList[[i]]}
                                 
              
              graphics::plot(x=df_qc[names(plotCols),"CellIndex"],y=y_i,pch=20,col=plotCols,
                             xlab="Cells ordered by lane",ylab=plotNames[i],cex.axis=1.2,cex.lab=1.2,
                             xlim=c(0,max(df_qc[names(plotCols),"CellIndex"])+5))
                if (thresholdList[[i]]) {
                  graphics::abline(h=thres_i,lty=lty,lwd=lwd)
                  graphics::text(x=max(df_qc[names(plotCols),"CellIndex"]), y=thres_i+0.1, thresholdList[[i]], col = "blue") 
                  
                }
                if (plotting == "pdf"){
                  
                  grDevices::pdf(paste0(outputPlots, "qc_plot_", plotList[[i]], ".pdf"), width = 7, height = 5)
                  
                  graphics::plot(x=df_qc[names(plotCols),"CellIndex"],y=y_i,pch=20,col=plotCols,
                       xlab="Cells ordered by lane",ylab=plotNames[i],cex.axis=1.2,cex.lab=1.2,
                       xlim=c(0,max(df_qc[names(plotCols),"CellIndex"])+5))
                  if (thresholdList[[i]]) {
                    graphics::abline(h=thres_i,lty=lty,lwd=lwd)
                    graphics::text(x=max(df_qc[names(plotCols),"CellIndex"]), y=thres_i+0.1, thresholdList[[i]], col = "blue") 
                    
                  }
                  dev.off()
                } else if (plotting == "tiff"){
                  
                  grDevices::tiff(paste0(outputPlots, "qc_plot_", plotList[[i]], ".tiff"),width = 7, height = 5, units = 'in', res = 300)

                  graphics::plot(x=df_qc[names(plotCols),"CellIndex"],y=y_i,pch=20,col=plotCols,
                                 xlab="Cells ordered by lane",ylab=plotNames[i],cex.axis=1.2,cex.lab=1.2,
                                 xlim=c(0,max(df_qc[names(plotCols),"CellIndex"])+5))
                  if (thresholdList[[i]]) {
                    graphics::abline(h=thres_i,lty=lty,lwd=lwd)
                    graphics::text(x=max(df_qc[names(plotCols),"CellIndex"]), y=thres_i+0.1, thresholdList[[i]], col = "blue") 
                    
                  }
                  dev.off()
                  
                  
                }
                
              }
              
            
            
            
            
            
            # Output cells passing/failingQC ------------------------------------------
            #geneExpressionMatrix <- counts[nuclearGeneIDs,!failQCBool]
            
            
            cat("DONE!\n\n")
            return(failQCBool)
            
              }
              )






