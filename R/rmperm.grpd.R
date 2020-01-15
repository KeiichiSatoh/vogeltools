#' Randomly Permute the Rows and Columns of an Input Matrix by group of actors.
#'
#' Given an input matrix (or stack thereof), \code{rmperm.grpd} performs a
#' (random) simultaneous row/column permutation of the input data BY EACH GROUP.
#' That means permutation does not occur among the different groups.
#'
#' This is the extention of the \code{\link[sna]{rmperm}} for the matrix
#' permtation. If no group was defined in \code{grp.ns}, it performs the normal
#' matrix permutation. If groups are defined, it perfom the matrix permutation
#' by each group seperately and then convert it into a single (or stacked)
#' matrix. This function helps for QAP test in which the off-diagonal matrix is
#' set to be \code{NA}, for example.
#'
#' @param mat a matrix, or stack thereof. For a stack of matrices, it should be either list format or array in which the first dimension indicates the layer.
#' @param grp.ns a vector indicates the number of actors in each group. The default is \code{NULL}.
#' @return a permuted matrix (or  matrices).
#'
#' @examples
#' library(sna)
#' # single group/single net(normal permutation)
#' m <- sna::rgraph(4)
#' rmperm.grpd(m)
#'
#' # single group/multiplex(normal permutation)
#' m <- rgraph(4,2)
#' m[1,,];m[2,,]
#' res <- rmperm.grpd(m)
#' res[1,,];res[2,,]
#'
#' # multi group/single net
#' m1 <- sna::rgraph(4)
#' m2 <- sna::rgraph(3)
#' m3 <- sna::rgraph(3)
#' spr.mat <- list2mat(list(m1,m2,m3))
#' spr.mat
#' rmperm.grpd(spr.mat, grp.ns = c(4,3,3))
#'
#' # multi group/multiplex
#' spr.mat2 <- spr.mat
#' spr.mat2[3,2]<-1;spr.mat2[1,3]<-1
#' sprmats <- list(spr.mat,spr.mat2)
#' res <- rmperm.grpd(sprmats, grp.ns = c(nrow(m1),nrow(m2),nrow(m3)))
#' res[1,,];res[2,,]
#'
#' @export

rmperm.grpd <- function(mat, grp.ns = NULL){
  # Change the condition by the group size
  if(is.null(grp.ns) == TRUE || length(grp.ns)==1){
    permed.mat <- sna::rmperm(mat)
    return(permed.mat)
  }
  # devide the graphs into groups
  partition <- cumsum(grp.ns)
  num.of.grp <- length(grp.ns)

  # Changing the condition of the network type
  mat <- sna::as.sociomatrix.sna(mat)

  ## single matrix
  if(length(dim(mat))==2){
    # make the permuted order for each group
    perm.order <- vector("list", num.of.grp)
    for(i in 1:num.of.grp){
      perm.order[[i]] <- sample(1:grp.ns[i])
    }
    # create a separate matrix
    grp.mat <- vector("list", num.of.grp)
    # First group
    grp.mat[[1]] <- mat[1:partition[1], 1:partition[1]]
    # Second group and so forth
    for(i in 2:num.of.grp){
      grp.mat[[i]] <- mat[(partition[i-1]+1):partition[i],
                        (partition[i-1]+1):partition[i]]
    }
    # Permute the matrix separately
    perm.grp.mat <- vector("list", num.of.grp)
    for(i in 1:num.of.grp){
      perm.grp.mat[[i]] <- grp.mat[[i]][perm.order[[i]],
                                        perm.order[[i]]]
    }
    # put back to the spread matrix
    permed.spr.mat <- list2mat(perm.grp.mat)
  }else{
    # make the permuted order for each group
    perm.order <- vector("list", num.of.grp)
    for(i in 1:num.of.grp){
      perm.order[[i]] <- sample(1:grp.ns[i])
    }
    # create a separate matrix
    grp.mat <- vector("list", num.of.grp)
    grp.mat.list <- vector("list", dim(mat)[1])
    for(k in 1:length(grp.mat.list)){
      grp.mat.list[[k]] <- grp.mat
    }
    # First group
    for(k in 1:length(grp.mat.list)){
      grp.mat.list[[k]][[1]] <- mat[k,
                                  1:partition[1], 1:partition[1]]
    }
    # Second group and so forth
    for(k in 1:length(grp.mat.list)){
      for(i in 2:num.of.grp){
        grp.mat.list[[k]][[i]] <- mat[k,
                                    (partition[i-1]+1):partition[i],
                                    (partition[i-1]+1):partition[i]]
      }
    }
    # Permute the matrix separately
    perm.grp.mat <- grp.mat.list
    for(k in 1:length(grp.mat.list)){
      for(i in 1:num.of.grp){
        perm.grp.mat[[k]][[i]] <- grp.mat.list[[k]][[i]][perm.order[[i]],
                                                         perm.order[[i]]]
      }
    }
    # put back to the spread matrix
    permed.spr.mat <- vector("list",length(grp.mat.list))
    for(k in 1:length(grp.mat.list)){
      permed.spr.mat[[k]] <- list2mat(perm.grp.mat[[k]])
    }
    # attach in the array form
    permed.spr.mat <- sna::as.sociomatrix.sna(permed.spr.mat)
  }
  return(permed.spr.mat)
}
