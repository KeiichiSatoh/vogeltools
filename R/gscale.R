#' Scaling and centering of matrix.
#'
#' \code{gscale} standardize the values in a matrix. This is a matrix version of
#' \code{scale} in the R package \code{base}.
#'
#' This function first vectrizes the values in a matrix and then passes to the
#' \code{scale} function in R package  \code{base}. The returned values are then
#' converted into the original matrix format. If \code{center} is \code{TRUE},
#' then centering is done by subtracting the means (omitting \code{NA}s) of the
#' matrix from the original values in the matrix. If \code{scale} is
#' \code{TRUE}, then scaling is done by dividing the (centered) values in the
#' matrix by their standard deviations if \code{center} is \code{TRUE}, and the
#' root mean square otherwise. If \code{scale} is \code{FALSE}, no scaling is done.
#' If \code{diag} is \code{TURE}, the diagonal values are included for the
#' transformation.
#'
#' @param mat a matrix to be transformed.
#' @param center a logical value indicating whether centering should be done.
#'   The default value is \code{TRUE}.
#' @param scale a logical value indicating whether scaling should be done. The
#'   default value is \code{TRUE}.
#' @param diag Whether the diagonal value is considered or not. The defalt value
#'   is \code{FALSE}(not considered). With this defalut setting, the diagonal
#'   values are set to be \code{NA}.
#' @return a scaled matrix defined by the scaling method.
#'
#' @examples
#' mat <- matrix(c(1,2,3,
#'                4,5,6,
#'                7,8,9),3,3,byrow = TRUE)
#' dimnames(mat) <- list(1:3, c("A","B","C"))
#'
#' gscale(mat)
#' gscale(mat, diag = TRUE)
#' gscale(mat, center = TRUE, scale = FALSE, diag = TRUE)
#' gscale(mat, center = FALSE, scale = TRUE, diag = TRUE)
#'
#' @export

gscale <- function(mat, center = TRUE, scale = TRUE, diag = FALSE){
  dim.mat <- dim(mat)
  if(diag == FALSE){
    diag(mat) <- NA
  }
  gmat <- matrix(scale(as.vector(mat),
                       center = center,
                       scale = scale),
                 dim.mat[1], dim.mat[2])
  dimnames(gmat) <- dimnames(mat)
  return(gmat)
}
