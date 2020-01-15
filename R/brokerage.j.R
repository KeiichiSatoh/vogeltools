#' Positional analysis of j and k in the Brokerage relationship i-j-k.
#'
#' \code{brokerage.j} and \code{brokerage.k} counts the number of times the
#' alter stands at the position in j and k, respectively, in the different
#' brokerage typology
#'
#' This function first differentiate the ties into homophilous and
#' heterophilious ties. It counts then how many times the alter stands in the j
#' position in the i-j-k relationship from the perspective of i. In the result
#' matrix the value of the cell M[i,j] indicates the number of times the alter j
#' stands in the j position. The count is standardized in accordance to the
#' method defined in \code{denom}. Brokerage analysis was first introduced in
#' Gould and Fernandez (1989). This function modifies the original concept in a
#' way that we construct the brokerage typology not based on the subgroup that
#' consists of actors engaging in the same issues, but the homophily and
#' heterophily ties.
#'
#' In the brokerage relationship where i and j as well as j and k are connected
#' AND i and k is not connected,
#'
#' Coordinator (wi) is the relationship where i->j and j->k is homophilous ties;
#'
#' Itinerant broker (wo) is the one where i->j and j->k is heterophilious ties;
#'
#' Representative (bio) is the one where i->j is a homophilous and j->k is a
#' heterophilious tie;
#'
#' Gatekeeper (boi) is the one where i->j is a heterophilious and j->k is a
#' homophilous tie.
#'
#' The fifth type, liaison, in the original literature (Gould and Fernandez
#' 1989) does not come up in this function due to the aforementioned
#' modification: the heterophiles tie between i->j and j->k already implies the
#' same type to itinerant broker(bo).
#'
#' The returned matrix M[i,j] is standardized if it is selected in
#' \code{"denom"}. \code{"n"} returns the value M[i,j]/n where n is number of
#' actors in the \code{mat}. \code{"density"} returns M[i,j]*density(mat).
#'
#' @param mat a matrix of a graph.
#' @param attr a vector of attributes.
#' @param method \code{"abs.diff"} or \code{"match"}. This is a intiger indicating the
#'   method how to calculate the difference of the \code{attr} between i and j.
#'   The default method \code{"abs.diff"} take the absolute
#'   abs(attr[i]-attr[j]). \code{"match"} code 1 if attr[i]==attr[j] and 0
#'   otherwise.
#' @param alpha a cut-off value to differentiate homophilous and heterophilious
#'   ties. If the difference of \code{attr} values between i and j exceeds the
#'   value of \code{alpha}, then the relationship between i and j is treated as
#'   heterophilious relationship. Otherwise homophilous.
#' @param denom \code{"none"}(default), \code{"n"} or \code{"density"}. A
#'   integer indicating how to standardize the counted number of times an actor
#'   stands at the j position.
#' @return a list of matrices in which each matrix shows the result of each
#'   typology.
#' @references Gould, Roger V. and Roberto M. Fernandez. (1989). "Structures of
#'   Mediation: A Formal Approach to Brokerage in Transaction Networks,"
#'   Sociological Methodology 19: 89-126.
#' @examples
#' library(sna)
#' mat <- sna::rgraph(5, tprob = 0.5, mode = "graph")
#' dimnames(mat) <- list(LETTERS[1:5],LETTERS[1:5])
#' attr <- c(1,1,1,0,0)
#' gplot(mat, displaylabels = TRUE, vertex.col = attr, usearrows = FALSE)
#'
#' brokerage.j(mat, attr, "match")
#' brokerage.j(mat, attr, "match", denom = "n")
#'
#' @export

brokerage.j <- function(
  mat = mat, attr = attr,
  method = "abs.diff", alpha = 0.5,
  denom = "none"){
  # calculate the difference of the attribute values
  if(method == "abs.diff"){
    diff.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
    dimnames(diff.mat) <- dimnames(mat)
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        diff.mat[i,j] <- abs(attr[i]-attr[j])
      }
    }
    # differentiate the matrix
    Homo <- mat
    Homo[diff.mat > alpha] <- 0
    Hetr <- mat
    Hetr[diff.mat <= alpha] <- 0
  }else if(method=="match"){
    match.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
    dimnames(match.mat) <- dimnames(mat)
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        match.mat[i,j] <- attr[i]==attr[j]
      }
    }
    # differentiate the matrix
    Homo <- mat
    Homo[match.mat==F] <- 0
    Hetr <- mat
    Hetr[match.mat==T] <- 0
  }

  # Calculate the brokerage
  # wi
  wi.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  wi <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(wi) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          wi.calc[i,j,k] <- NA
        }else{
          wi.calc[i,j,k] <- (Homo[i,j]==1 & Homo[j,k]==1 & Homo[i,k]==0)
        }
        wi[i,j] <- sum(wi.calc[i,j,], na.rm = T)
      }
    }
  }

  # wo
  wo.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  wo <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(wo) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          wo.calc[i,j,k] <- NA
        }else{
          wo.calc[i,j,k] <- (Hetr[i,j]==1 & Hetr[j,k]==1 & Homo[i,k]==0)
        }
        wo[i,j] <- sum(wo.calc[i,j,], na.rm = T)
      }
    }
  }

  # boi
  boi.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  boi <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(boi) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          boi.calc[i,j,k] <- NA
        }else{
          boi.calc[i,j,k] <- (Hetr[i,j]==1 & Homo[j,k]==1 & Hetr[i,k]==0)
        }
        boi[i,j] <- sum(boi.calc[i,j,], na.rm = T)
      }
    }
  }

  # bio
  bio.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  bio <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(bio) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          bio.calc[i,j,k] <- NA
        }else{
          bio.calc[i,j,k] <- (Homo[i,j]==1 & Hetr[j,k]==1 & Hetr[i,k]==0)
        }
        bio[i,j] <- sum(bio.calc[i,j,], na.rm = T)
      }
    }
  }

  # return the result
  if(denom == "none"){
    res <- list(wi = wi,
                wo = wo,
                boi = boi,
                bio = bio)
    return(res)
  }else if(denom == "n"){
    n <- nrow(mat)
    res <- list(wi = wi/n,
                wo = wo/n,
                boi = boi/n,
                bio = bio/n)
    return(res)
  }else if(denom == "density"){
    D <- sna::gden(mat)
    res <- list(wi = wi*D,
                wo = wo*D,
                boi = boi*D,
                bio = bio*D)
    return(res)
  }
}

#' @rdname brokerage.j
brokerage.k <- function(
  mat = mat, attr=attr,
  method="abs.diff", alpha = 0.5,
  denom = "none"){
  if(method == "abs.diff"){
    diff.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
    dimnames(diff.mat) <- dimnames(mat)
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        diff.mat[i,j] <- abs(attr[i]-attr[j])
      }
    }
    # differentiate the matrix
    Homo <- mat
    Homo[diff.mat > alpha] <- 0
    Hetr <- mat
    Hetr[diff.mat <= alpha] <- 0
  }else if(method=="match"){
    match.mat <- matrix(NA, nrow(mat), ncol(mat), byrow = T)
    dimnames(match.mat) <- dimnames(mat)
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        match.mat[i,j] <- attr[i]==attr[j]
      }
    }
    # differentiate the matrix
    Homo <- mat
    Homo[match.mat==F] <- 0
    Hetr <- mat
    Hetr[match.mat==T] <- 0
  }

  # Calculate the brokerage
  # wi
  wi.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  wi <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(wi) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          wi.calc[i,j,k] <- NA
        }else{
          wi.calc[i,j,k] <- (Homo[i,j]==1 & Homo[j,k]==1 & Homo[i,k]==0)
        }
        wi[i,k] <- sum(wi.calc[i,,k], na.rm = T)
      }
    }
  }

  # wo
  wo.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  wo <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(wo) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          wo.calc[i,j,k] <- NA
        }else{
          wo.calc[i,j,k] <- (Hetr[i,j]==1 & Hetr[j,k]==1 & Homo[i,k]==0)
        }
        wo[i,k] <- sum(wo.calc[i,,k], na.rm = T)
      }
    }
  }

  # boi
  boi.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  boi <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(boi) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          boi.calc[i,j,k] <- NA
        }else{
          boi.calc[i,j,k] <- (Hetr[i,j]==1 & Homo[j,k]==1 & Hetr[i,k]==0)
        }
        boi[i,k] <- sum(boi.calc[i,,k], na.rm = T)
      }
    }
  }

  # bio
  bio.calc <- array(NA, c(nrow(mat), nrow(mat), nrow(mat)))
  bio <- matrix(NA, nrow(mat), nrow(mat), byrow = T)
  dimnames(bio) <- dimnames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat)){
      for(k in 1:nrow(mat)){
        if(i==j | j==k | k==i){
          bio.calc[i,j,k] <- NA
        }else{
          bio.calc[i,j,k] <- (Homo[i,j]==1 & Hetr[j,k]==1 & Hetr[i,k]==0)
        }
        bio[i,k] <- sum(bio.calc[i,,k], na.rm = T)
      }
    }
  }

  # return the result
  if(denom == "none"){
    res <- list(wi = wi,
                wo = wo,
                boi = boi,
                bio = bio)
    return(res)
  }else if(denom == "n"){
    n <- nrow(mat)
    res <- list(wi = wi/n,
                wo = wo/n,
                boi = boi/n,
                bio = bio/n)
    return(res)
  }else if(denom == "density"){
    D <- sna::gden(mat)
    res <- list(wi = wi*D,
                wo = wo*D,
                boi = boi*D,
                bio = bio*D)
    return(res)
  }
}

