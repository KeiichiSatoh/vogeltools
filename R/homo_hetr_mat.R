#' Differentiate homophilous and heterophilious graph.
#'
#' \code{homo_hetr_mat} takes a matrix and returns a list of matrices in which
#' homophilous and heterophilious ties are differentiated.
#'
#' @param mat a matrix of a graph.
#' @param attr a vector of attributes.
#' @param alpha a cut-off value to differentiate homophilous and heterophilious
#'   ties. If the difference of \code{attr} values between i and j exceeds the
#'   value of \code{alpha}, then the relationship between i and j is treated as
#'   heterophilious relationship. Otherwise homophilous.
#' @return a list of matrices: a homopilious matrix and a heterophilious matrix.
#' @examples
#' library(sna)
#' mat <- sna::rgraph(5)
#' dimnames(mat) <- list(LETTERS[1:5], LETTERS[1:5])
#' attr <- c(0,0,0,1,1)
#'
#' homo_hetr_mat(mat, attr)
#'
#' @export

homo_hetr_mat <- function(mat, attr, alpha = 0.5){
  diff.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
  dimnames(diff.mat) <- dimnames(mat)
  for(i in 1:nrow(mat)){for(j in 1:nrow(mat)){
    diff.mat[i,j] <- abs(attr[i] - attr[j])
  }}

  Homo <- mat
  Homo[diff.mat>alpha] <- 0
  Hetr <- mat
  Hetr[diff.mat<=alpha] <- 0

  res <- list(Homo=Homo,
              Hetr=Hetr)
  return(res)
}
