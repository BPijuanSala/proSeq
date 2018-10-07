########################################################################################
## Title: function_quantileNorm.R
## Author: Blanca Pijuan-Sala
## Description: Method to perform quantile normalisation.
## Date: 24 December 2016
## **proSeq package**
########################################################################################



#' @title Quantile/Rank normalisation
#' @description It normalises based on ranks.
#' @param data Dataframe or matrix containing the raw counts (genes x cells).
#' @param ties.rank It specifies how ties are treated. See \code{\link[base]{rank}} for
#' further details. Default: min.
#' @examples 
#' 
#' #test the function
#' data <- data.frame(one=c(5,2,3,4),
#'                    two=c(4,1,4,2),
#'                    three=c(3,4,6,8)
#' )
#' rownames(data) <- toupper(letters[1:4])
#' quantileNorm(data)

#' @references https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#' @author Blanca Pijuan-Sala.
#' @export
#' @rdname quantileNorm

quantileNorm <- function(data,ties.rank="min"){
  data_rank <- apply(data,2,rank,ties.method=ties.rank)
  data_sorted <- data.frame(apply(data, 2, sort))
  data_mean <- apply(data_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  data_final <- apply(data_rank, 2, index_to_mean, my_mean=data_mean)
  rownames(data_final) <- rownames(data)
  return(data_final)
}
