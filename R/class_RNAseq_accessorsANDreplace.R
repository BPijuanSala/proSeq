########################################################################################
## Title: class_RNAseq_accessorsANDreplace.R
## Author: Blanca Pijuan-Sala
## Description: Definition of RNAseq class and their accessor and replacement methods.
## Date: 20 December 2016
## **proSeq package**
########################################################################################



##############################################################
## Class RNAseq
##############################################################

#' @title Class RNAseq
#' @description Class \code{RNAseq} defines an RNAseq dataset.
#' @name RNAseq-class
#' @rdname RNAseq-class
#' @slot countsRaw matrix of raw counts (genes x cells). Recommended: Use gene IDs 
#' instead of gene Names.
#' @slot CountsFeat matrix (features x cells) indicating the number of 
#' no feature reads (row 1), ambiguous reads (row 2), low quality reads (row 3),
#' unaligned reads (row 4) and others. The cellNames should be the same as in countsRaw.
#' @slot metadata matrix of metadata of all cells and genes (cells x features).
#' @slot SpikeIn If not NULL, contains logical vector with genes as names, 
#' where TRUE indicates that it's spike-in and FALSE indicates that it's a biological gene.
#' @slot countsNorm matrix of normalised counts (genes x cells).Recommended: 
#' Use gene IDs instead of gene Names.
#' @slot CellsSizeFac Contains size factors obtained in the normalisation.
#' @slot cellsFailQC List of cells that have failed QC.
#' @slot genesHVG Vector of highly variable genes.
#' @slot pca If not NULL, contains PCA results.
#' @slot tsne If not NULL, contains tSNE results.
#' @slot difmap If not NULL, contains diffusion map results.
#' @slot clustering If not NULL, contains information about potential clustering.
#' @exportClass RNAseq

setClass("RNAseq",slots = c(
  countsRaw="matrix",
  CountsFeat="matrix",
  metadata="matrix",
  SpikeIn="logical",
  countsNorm="matrix",
  CellsSizeFac="list",
  #it must include (1) sizeFactors for biological genes and spike-ins.
  cellsFailQC="vector",#this includes cells that pass QC.
  genesHVG="vector",
  pca="ANY",
  #output of prcomp or printcomp
  tsne="list",
  #output of Rtsne
  difmap="ANY",
  #output of diffusion map
  clustering="ANY"
  #Apart from the clustering out, the slot must contain at least: (1) vector of cluster numbers and named by cells
  # and (2) type of clustering (simple vector with name), (3) parameters: vector with values 
  ## named by the parameter names.
  ),
  validity = function(object) {
    isValid <- TRUE
    if (((dim(object@countsNorm)[1] == 0) & (dim(object@metadata)[1]==0))==FALSE) {
      if (length(which((colnames(object@countsNorm) %in% rownames(object@metadata))==FALSE)) > 0) {
        isValid <- FALSE
        cat("Not all the cells in countsNorm have metadata. Please include metadata for all of them. You can 
            fill the metadata as NA if necessary.")
        #cat("The number of cells in the metadata (rows) is not the same as the number
        #    of cells in the counts matrix (columns). Please introduce all the cells for these 
        #    two slots.\n")
        
      }
      if (ncol(object@countsNorm) > nrow(object@metadata)) {
        isValid <- FALSE
        cat("Not all the cells in countsNorm have metadata.Please add metadata for them. you can include NA for the cells you lack information.")
        #cat("The number of cells in the metadata (rows) is not the same as the number
        #    of cells in the counts matrix (columns). Please introduce all the cells for these 
        #    two slots.\n")
      
      }
      
      
    }
    
    if (((dim(object@countsRaw)[1] == 0) & (dim(object@metadata)[1]==0))==FALSE) {
      if (length(which((colnames(object@countsRaw) %in% rownames(object@metadata))==FALSE)) > 0) {
        isValid <- FALSE
        cat("Not all the cells in countsRaw have metadata.Please add metadata for them. you can include NA for the cells you lack information.")
        #cat("The number of cells in the metadata (rows) is not the same as the number
        #    of cells in the counts matrix (columns). Please introduce all the cells for these 
        #    two slots.\n")
        
      }
  
    return(isValid)
    }
  }
  )




####################################################################
## Accessor and replacement methods
####################################################################


#' @title Show details of \code{RNAseq} class
#' @description Show details of \code{RNAseq} class.
#' @param object \code{RNAseq} object.
#' @return It gives general information on the RNAseq class.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname show_RNAseq-method
setGeneric("show", function(object) standardGeneric("show"))

#' @rdname show_RNAseq-method
#' @export
setMethod("show",
          "RNAseq",
          function(object) {
            cat("Object of class",class(object),".\n")
            cat("It contains space for:\n")
            cat("  countsRaw (",class(object@countsRaw),"): matrix of raw counts (genes x cells).\n")
            cat("  CountsFeat (",class(object@CountsFeat),"): matrix of features x cells indicating the number of 
                  feature reads (column 1), ambiguous reads (column 2), low quality reads (column 3), 
                unaligned reads (column 4). The cellNames should be the same as in countsRaw. \n")
            cat("  metadata (",class(object@metadata),"): matrix of metadata of all cells and genes (cells x features).\n")
            cat("  SpikeIn (", class(object@SpikeIn),"): If not NULL, contains logical vector with genes as names, 
            where TRUE indicates that it's spike-in and FALSE indicates that it's a biological gene. \n")
            cat("  countsNorm (",class(object@countsNorm),"): matrix of normalised counts (genes x cells). \n")
            cat("  CellsSizeFac (",class(object@countsNorm),"): Contains size factors obtained in the normalisation.\n")
            cat("  cellsFailQC (",class(object@cellsFailQC),"): List of cells that have passed QC. \n")
            cat("  genesHVG (",class(object@genesHVG),"): Vector of highly variable genes. \n")
            cat("  pca (",class(object@pca),"): If not NULL, contains PCA results. \n")
            cat("  tsne (",class(object@tsne),"):If not NULL, contains tSNE results. \n")
            cat("  difmap (",class(object@difmap),"): If not NULL, contains diffusion map results. \n")
            cat("  clustering (",class(object@clustering),"): If not NULL, contains information about potential clustering. \n")
            
          })


#' @title Access countsRaw data
#' @description Accesses \code{countsRaw} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return countsRaw
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname countsRaw_access
setGeneric("countsRaw", function(object) standardGeneric("countsRaw"))

#' @rdname countsRaw_access
#' @export
setMethod("countsRaw", "RNAseq", function(object) object@countsRaw)


#' @title Modify countsRaw data
#' @description Modifies \code{countsRaw} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Matrix or dataframe containing new countsRaw data for 
#' \code{RNAseq} object.
#' @return It does not display any value but modifies countsRaw.
#' @author Blanca Pijuan Sala
#' @export
#' @rdname countsRaw_replace
setGeneric("countsRaw<-", function(object,value) standardGeneric("countsRaw<-"))

#' @rdname countsRaw_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("countsRaw",signature(object="RNAseq",value="matrix"),
                 function(object,value){
                   object@countsRaw <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)



#' @title Access metadata
#' @description Accesses \code{metadata} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return metadata
#' @author Blanca Pijuan Sala
#' @export
#' @rdname metadata_access
setGeneric("metadata", function(object) standardGeneric("metadata"))

#' @rdname metadata_access
#' @export
setMethod("metadata", "RNAseq", function(object) object@metadata)



#' @title Modify metadata slot
#' @description Modifies \code{metadata} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Matrix or dataframe containing new metadata for \code{RNAseq} object.
#' @return It does not display any value but modifies metadata.
#' @author Blanca Pijuan Sala
#' @export
#' @rdname metadata_replace
setGeneric("metadata<-", function(object,value) standardGeneric("metadata<-"))

#' @rdname metadata_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("metadata",signature(object="RNAseq",value="matrix"),
                 function(object,value){
                   object@metadata <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)



#' @title Access countsNorm data
#' @description Accesses \code{countsNorm} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return countsNorm.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname countsNorm_access
setGeneric("countsNorm", function(object) standardGeneric("countsNorm"))

#' @rdname countsNorm_access
#' @export
setMethod("countsNorm", "RNAseq", function(object) object@countsNorm)



#' @title Modify countsNorm data
#' @description Modifies \code{countsNorm} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Matrix or dataframe containing new countsNorm data
#' for \code{RNAseq} object.
#' @return It does not display any value but modifies countsNorm.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname countsNorm_replace
setGeneric("countsNorm<-", function(object,value) standardGeneric("countsNorm<-"))

#' @rdname countsNorm_replace
#' @export
#' @importFrom methods validObject

setReplaceMethod("countsNorm",signature(object="RNAseq",value="matrix"),
                 function(object,value){
                   object@countsNorm <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)



#' @title Access CountsFeat data
#' @description Accesses \code{CountsFeat} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return CountsFeat.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname CountsFeat_access
setGeneric("CountsFeat", function(object) standardGeneric("CountsFeat"))

#' @rdname CountsFeat_access
#' @export
setMethod("CountsFeat", "RNAseq", function(object) object@CountsFeat)



#' @title Modify CountsFeat data
#' @description Modifies \code{CountsFeat} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Matrix or dataframe containing new CountsFeat data
#' for \code{RNAseq} object.
#' @return It does not display any value but modifies CountsFeat.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname CountsFeat_replace
setGeneric("CountsFeat<-", function(object,value) standardGeneric("CountsFeat<-"))

#' @rdname CountsFeat_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("CountsFeat",signature(object="RNAseq",value="matrix"),
                 function(object,value){
                   object@CountsFeat <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)





#' @title Access CellsSizeFac data
#' @description Accesses \code{CellsSizeFac} slot in \code{RNAseq} class object
#' @param object \code{RNAseq} object.
#' @return CellsSizeFac
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname CellsSizeFac_access
setGeneric("CellsSizeFac", function(object) standardGeneric("CellsSizeFac"))

#' @rdname CellsSizeFac_access
#' @export
setMethod("CellsSizeFac", "RNAseq", function(object) object@CellsSizeFac)


#' @title Modify CellsSizeFac data
#' @description Modifies \code{CellsSizeFac} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new CellsSizeFac for \code{RNAseq} object.
#' @return It does not display any value but modifies CellsSizeFac.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname CellsSizeFac_replace
setGeneric("CellsSizeFac<-", function(object,value) standardGeneric("CellsSizeFac<-"))

#' @rdname CellsSizeFac_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("CellsSizeFac",signature(object="RNAseq",value="list"),
                 function(object,value){
                   object@CellsSizeFac <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)

#' @title Access cellsFailQC data
#' @description Accesses \code{cellsFailQC} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return cellsFailQC.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname cellsFailQC_access
setGeneric("cellsFailQC", function(object) standardGeneric("cellsFailQC"))

#' @rdname cellsFailQC_access
#' @export
setMethod("cellsFailQC", "RNAseq", function(object) object@cellsFailQC)


#' @title Modify cellsFailQC data
#' @description Modifies \code{cellsFailQC} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new cellsFailQC data for \code{RNAseq} object.
#' @return It does not display any value but modifies cellsFailQC.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname cellsFailQC_replace
setGeneric("cellsFailQC<-", function(object,value) standardGeneric("cellsFailQC<-"))

#' @rdname cellsFailQC_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("cellsFailQC",signature(object="RNAseq",value="vector"),
                 function(object,value){
                   object@cellsFailQC <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)



#' @title Access genesHVG data
#' @description Accesses \code{genesHVG} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return genesHVG.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname genesHVG_access
setGeneric("genesHVG", function(object) standardGeneric("genesHVG"))

#' @rdname genesHVG_access
#' @export
setMethod("genesHVG", "RNAseq", function(object) object@genesHVG)


#' @title Modify genesHVG data
#' @description Modifies \code{genesHVG} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new genesHVG data for \code{RNAseq} object.
#' @return It does not display any value but modifies genesHVG.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname genesHVG_replace
setGeneric("genesHVG<-", function(object,value) standardGeneric("genesHVG<-"))

#' @rdname genesHVG_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("genesHVG",signature(object="RNAseq",value="vector"),
                 function(object,value){
                   object@genesHVG <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)




#' @title Access pca data
#' @description Accesses \code{pca} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return pca.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname pca_access
setGeneric("pca", function(object) standardGeneric("pca"))

#' @rdname pca_access
#' @export
setMethod("pca", "RNAseq", function(object) object@pca)


#' @title Modify pca data
#' @description Modifies \code{pca} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new pca data for \code{RNAseq} object.
#' @return It does not display any value but modifies pca.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname pca_replace
setGeneric("pca<-", function(object,value) standardGeneric("pca<-"))

#' @rdname pca_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("pca",signature(object="RNAseq",value="ANY"),
                 function(object,value){
                   object@pca <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)


#' @title Access tsne data
#' @description Accesses \code{tsne} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return tsne.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname tsne_access
setGeneric("tsne", function(object) standardGeneric("tsne"))

#' @rdname tsne_access
#' @export
setMethod("tsne", "RNAseq", function(object) object@tsne)


#' @title Modify tsne data
#' @description Modifies \code{tsne} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new tsne data for \code{RNAseq} object.
#' @return It does not display any value but modifies tsne.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname tsne_replace
setGeneric("tsne<-", function(object,value) standardGeneric("tsne<-"))

#' @rdname tsne_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("tsne",signature(object="RNAseq",value="ANY"),
                 function(object,value){
                   object@tsne <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)


#' @title Access difmap data
#' @description Accesses \code{difmap} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return difmap.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname difmap_access
setGeneric("difmap", function(object) standardGeneric("difmap"))

#' @rdname difmap_access
#' @export
setMethod("difmap", "RNAseq", function(object) object@difmap)



#' @title Modify difmap data
#' @description Modifies \code{difmap} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new difmap data for \code{RNAseq} object.
#' @return It does not display any value but modifies difmap.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname difmap_replace
setGeneric("difmap<-", function(object,value) standardGeneric("difmap<-"))

#' @rdname difmap_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("difmap",signature(object="RNAseq",value="ANY"),
                 function(object,value){
                   object@difmap <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)


#' @title Access clustering data
#' @description Accesses \code{clustering} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @return clustering.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname clustering_access
setGeneric("clustering", function(object) standardGeneric("clustering"))

#' @rdname clustering_access
#' @export
setMethod("clustering", "RNAseq", function(object) object@clustering)


#' @title Modify clustering data
#' @description Modifies \code{clustering} slot in \code{RNAseq} class object.
#' @param object \code{RNAseq} object.
#' @param value Object containing new clustering data for \code{RNAseq} object.
#' @return It does not display any value but modifies clustering.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname clustering_replace
setGeneric("clustering<-", function(object,value) standardGeneric("clustering<-"))

#' @rdname clustering_replace
#' @export
#' @importFrom methods validObject
setReplaceMethod("clustering",signature(object="RNAseq",value="ANY"),
                 function(object,value){
                   object@clustering <- value
                   if (validObject(object)){
                     return(object)} else {
                       cat("input not valid.")
                     }
                 }
)



#run_if = FALSE
##Some checks:
#if (run_if == TRUE){
#  #trying:
#  num_cells = 20
#  counts_try= matrix(rexp(200, rate=.1), ncol=num_cells)
#  colnames(counts_try)=paste0("cell_",seq(1:ncol(counts_try)))
#  metadata_try= matrix(rexp(200, rate=.1), nrow=num_cells)
  #rownames(metadata_try)[1] <- "jifrjdksf"
##  rownames(metadata_try)=paste0("cell_",seq(1:nrow(metadata_try)))

 # data_try<- new("RNAseq", countsRaw=counts_try, metadata=metadata_try)
#  pca_try=prcomp(counts_try)

#  pca(data_try) <- pca_try

  
  #tsne_try=Rtsne(counts_try,perplexity=1)
  #data_try@tsne <- tsne_try

 # dif_map_try=DiffusionMap(counts_try)
#  data_try@difmap <- dif_map_try
#}
#setClass("A", representation(x="numeric"))
#setClass("B", representation(y="numeric"), contains="A")

