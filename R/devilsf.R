#' Create a devil shift and angel shift matrix.
#'
#' \code{devilsf} returns the devil shift and angel shift matrix
#'
#' Two devil shift and angel shift matrix is created from the two scale:(1) the
#' degree that the ego perceives the others more powerful than they actually are
#' (overestimation of powerfulness) and (2) the degree of opponents (devilness)
#' as well as proponents (angelness) from the perspective of respondents (egos).
#' Several options are available that need to be adjusted based on the
#' assumption of the researcher.
#'
#' @param mat the matrix showing the influence checked by the row actors to
#'   column actors.
#' @param pol the scale indicating the policy preference.
#' @param score.std string indicating the method of standardizing the policy
#'   prefrence. score.std is "none" by default.
#' @param exp.method string indicating the method for calculating the degree of
#'   overestimation.
#' @param eval.method string indicating the method to evaluate the
#'   overestimation.
#' @param alpha cut-off value when binalize the devil shift and angel shift
#'   matrix. alpha is NULL by default.
#' @param score.norm logical; should the engel- and devil shift score normalized
#'   by the total number of influence checked in the matrix. score.norm is FALSE
#'   by default.
#' @examples
#' mat <- matrix(c(0,0,0,1,1,1,
#'                 1,1,1,0,0,0,
#'                 1,0,0,1,0,0,
#'                 1,1,1,0,0,0,
#'                 0,0,0,1,1,1,
#'                 1,0,0,1,0,0),6,6,byrow = TRUE)
#' namelabel <- c(paste(rep("A",3), 1:3, sep = ""),
#'                paste(rep("B",3), 1:3, sep = ""))
#' dimnames(mat) <- list(namelabel, namelabel)
#'
#' pol <- c(1,1,0.6,0.4,0,0)
#' names(pol) <- namelabel
#' res <- devilsf(mat, pol,
#'               score.std = "none",
#'               exp.method ="both",
#'               eval.method = "realized",
#'               alpha = NULL,
#'               score.norm = FALSE)
#'
#' res$original.mat
#' res$exp.mat
#' res$eval.exp.mat
#' res$pol.dis
#' res$pol.agg
#' res$devil.mat
#' res$angel.mat
#'
#' @export
devilsf <- function(
  mat = mat,
  pol = pol,
  score.std = "none",
  exp.method = "both",
  eval.method = "realized",
  alpha = NULL,
  score.norm = FALSE){

  # transform the policy score ------------------------
  if(score.std=="none"){
    pol <- pol
  }else if(score.std=="zero.one"){
    pol.zero.one <- (pol-min(pol))/(max(pol)-min(pol))
    pol <- pol.zero.one
  }else if(score.std=="z.std"){
    pol.z <- scale(pol, center = T, scale = T)
    pol.z.std <- (pol.z-min(pol.z))/(max(pol.z)-min(pol.z))
    pol <- pol.z.std
  }else if(score.std=="non.param"){
    pol.rank <- rank(pol, ties.method = "average")
    pol.rank.std <- (pol.rank-min(pol.rank))/(max(pol.rank)-min(pol.rank))
    pol <- pol.rank.std
  }else{
    stop("no score.std method available")
  }

  # null matrix ----------------------------------------
  exp.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
  dimnames(exp.mat) <- dimnames(mat)
  mat.rowsum <- rowSums(mat)
  mat.colsum <- colSums(mat)
  mat.allsum <- sum(mat.rowsum)

  # expected value -------------------------------------
  if(exp.method == "both"){
    for(i in 1:nrow(mat)){
      for(j in 1:ncol(mat)){
        exp.mat[i,j] <- (mat.rowsum[i] * mat.colsum[j])/mat.allsum
      }
    }
  }else if(exp.method == "column"){
    for(i in 1:nrow(mat)){
      exp.mat[i,] <- mat.colsum/nrow(mat)
    }
  }else if(exp.method == "row"){
    for(j in 1:ncol(mat)){
      exp.mat[,j] <- mat.rowsum/nrow(mat)
    }
  }else if(exp.method == "none"){
    exp.mat <- mat
  }else{
    stop("no method for exp.method available")
  }

  # replace the cell whose expected value exceeds 1 to 1
  exp.mat[exp.mat>=1]  <- 1
  exp.mat[exp.mat<=-1] <- -1

  # evaluation matrix
  if(exp.method == "none"){
    eval.exp.mat <- mat
  }else if(eval.method == "realized"){
    subtract.mat <- mat - exp.mat
    eval.exp.mat <- mat * subtract.mat
  }else if(eval.method == "all"){
    subtract.mat <- mat - exp.mat
    eval.exp.mat <- subtract.mat
  }else{
    stop("no exp.method method available")}


  # creating the policy preference matrix ---------------------
  ## policy disagreement matrix
  pol.dis <- matrix(NA, length(pol), length(pol), byrow = T)
  dimnames(pol.dis) <- dimnames(mat)
  for(i in 1:nrow(pol.dis)){for(j in 1:ncol(pol.dis)){
    pol.dis[i,j] <- abs(pol[i]-pol[j])
  }}
  pol.agg <- max(pol.dis)-pol.dis

  # Differentiate whether we will use the proportional of agreement or binary
  if(is.null(alpha) == FALSE){
    pol.dis[pol.dis>=alpha] <- 1
    pol.dis[pol.dis<alpha]  <- 0
    pol.agg[pol.agg>=alpha] <- 1
    pol.agg[pol.agg<alpha]  <- 0
  }
  # devil shift matrix
  devil.mat <- eval.exp.mat * pol.dis
  # angel shift matrix
  angel.mat <- eval.exp.mat * pol.agg

  # devil and angel shift score(whole)
  angel.score <- sum(angel.mat)
  devil.score <- sum(devil.mat)

  # devil and angel shift score(actor)
  angel.out <- rowSums(angel.mat)
  angel.in  <- colSums(angel.mat)
  devil.out <- rowSums(devil.mat)
  devil.in  <- colSums(devil.mat)

  # score normalization
  if(score.norm == T){
    n.infl <- sum(mat)
    angel.score <- angel.score/n.infl
    devil.score <- devil.score/n.infl
    angel.out <- angel.out/rowSums(mat)
    angel.in <- angel.in/colSums(mat)
    devil.out <- devil.out/rowSums(mat)
    devil.in <- devil.in/colSums(mat)
  }
  # result
  res <- list(
    original.mat = mat,
    angel.score = angel.score,
    devil.score = devil.score,
    angel.out = angel.out,
    angel.in = angel.in,
    devil.out = devil.out,
    devil.in = devil.in,
    exp.mat = exp.mat,
    eval.exp.mat = eval.exp.mat,
    devil.mat = devil.mat,
    angel.mat = angel.mat,
    pol.dis = pol.dis,
    pol.agg = pol.agg
  )
  # return
  print(c(Angel=res$angel.score, Devil=res$devil.score))
  invisible(res)
}
