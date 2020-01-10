#' Convert List of matrices to a matrix.
#'
#' \code{list2mat} converts a list of matrices to a matrix in which each orginal
#' matrix locates in the diagonal and off-diagonal cells are set to be NA (by
#' default).
#'
#' This function convertes a original list of matrices to a new matrix in which
#' each orginal matrix locates in the diagonal and of diagonal cells are set to
#' be NA (the filler can be set by changing \code{fill}. This conversion is
#' useful, for example, when running an QAP regression with the intention to
#' find out the coefficients across networks.
#'
#' @param mats the matrices stored in a \code{list} format. \code{array} is also
#'   possible to supply as input, if the first dimension of the array indicats the
#'   each network, and the second and third are the row and column of the matrix
#'   (this is the default format of package \code{sna}).
#' @param fill the value of the off-diagonal cells. The defalut value is
#'   \code{NA}.
#' @return a matrix converted from the list of matrices.
#'
#' @examples
#' ## library(sna)
#' mat1 <- rgraph(3)
#' dimnames(mat1) <- list(paste(rep("A", 3), 1:3, sep = ""),
#'                       paste(rep("A", 3), 1:3, sep = ""))
#' mat2 <- rgraph(4)
#' dimnames(mat2) <- list(paste(rep("B", 4), 1:4, sep = ""),
#'                       paste(rep("B", 4), 1:4, sep = ""))
#' mat3 <- rgraph(2)
#' dimnames(mat3) <- list(paste(rep("C", 2), 1:2, sep = ""),
#'                       paste(rep("C", 2), 1:2, sep = ""))
#' mats <- list(mat1, mat2, mat3)
#'
#' list2mat(mats, fill = NA)

list2mat <- function(mats = mat.list, fill = NA){
  # convert array to matrix, if the original data is that format
  if(is.array(mats)==T){
    mats2 <- list()
    for(i in 1:dim(mats)[1]){
      mats2[[i]] <- mats[i,,]
    }
    mats <- mats2
  }

  # number of matrix
  m <- length(mats)

  # attach the dimension of each matrix
  dimlist <- matrix(NA, m, 2, byrow = T)
  rownames(dimlist) <- 1:m
  colnames(dimlist) <- c("row", "column")
  for(k in 1:m){
    dimlist[k,] <- dim(mats[[k]])
  }

  # cumurative dimlist
  rownum.cum <- cumsum(dimlist[,"row"])
  colnum.cum <- cumsum(dimlist[,"column"])

  # create the namelist
  rownamelist <- unlist(lapply(mats, rownames))
  colnamelist <- unlist(lapply(mats, colnames))

  # Create a new matrix
  newmat <- matrix(fill, sum(dimlist[,"row"]), sum(dimlist[,"column"]), byrow = T)
  dimnames(newmat) <- list(rownamelist, colnamelist)
  ## first list
  newmat[1:rownum.cum[1], 1:colnum.cum[1]] <- mats[[1]]
  ## second list and more
  for(k in 2:m){
    newmat[(rownum.cum[k-1]+1):rownum.cum[k],
           (colnum.cum[k-1]+1):colnum.cum[k]] <- mats[[k]]
  }
  return(newmat)
}



