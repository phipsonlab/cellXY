#' Normalise a counts matrix to the median library size
#'
#' @param x a \code{DGEList} object or matrix of counts.
#' @param log logical, indicates whether the output should be on the log2 scale
#' or counts scale. Default is FALSE.
#' @param prior.count The prior count to add if the data is log2 normalised.
#' Default is a small count of 0.5.
#' @param lib.size a vector of library sizes to be used during the normalisation
#' step. Default is NULL and will be computed from the counts matrix.
#'
#' @return a matrix of normalised counts
#' @importFrom stats median
#' @author Belinda Phipson
#'
lognormCounts <- function(x, log = TRUE, prior.count = 0.5, lib.size)
  # Function to log normalise a counts matrix to median library size
  # Input is counts matrix
  # Genes as rows, cells as columns
  # Belinda Phipson
  # 26 February 2020
{
  x <- as.matrix(x)
  #lib.size <- colSums(x)
  M <- median(lib.size)
  if(log){
    prior.count.scaled <- lib.size/mean(lib.size)*prior.count
    lib.size <- lib.size + 2*prior.count.scaled
    return(log2(t((t(x)+prior.count.scaled)/lib.size*M)))
  }else{
    return(t(t(x)/lib.size*M))
  }
}

