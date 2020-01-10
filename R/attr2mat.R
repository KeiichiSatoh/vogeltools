#' Convert attribute vector to a matrix
#'
#' \code{attr2mat} converts a vector indicating attribute to a matrix.
#'
#' This function convertes a vector of matrices \eqn{attr} to a new matrix
#' \eqn{mat[i,j]} based on the definition set by \code{method}:
#'
#' \code{"match"}(matching the value): \eqn{Matrix[i,j] = 1} if \eqn{attr[i] = attr[j]}; 0 otherwise.
#'
#' \code{"match.one"}(maching the value 1): \eqn{Matrix[i,j] = 1} if \eqn{attr[i] = attr[j] = 1}; 0 otherwise.
#'
#' \code{"diff"}(difference): \eqn{Matrix[i,j] = attr[i] - [j]}.
#'
#' \code{"abs.diff"}(absolute difference): \eqn{Matrix[i,j] = |attr[i] - attr[j]|}.
#'
#' \code{"sqd.diff"}(squared difference): \eqn{Matrix[i,j] = (attr[i] - attr[j]=1)^2}.
#'
#' \code{"sum"}(sum): \eqn{Matrix[i,j] = attr[i] + attr[j] = 1}.
#'
#' \code{"product"}(product): \eqn{Matrix[i,j] = attr[i] * attr[j]}.
#'
#' \code{"identity"}(identity coefficient): \eqn{Matrix[i,j] = 2 * (attr[i] * attr[j])/(attr[i]^2 + attr[j]^2}.
#'
#' \code{"row"}(duplicating row actors' attributes: sender effect): \eqn{Matrix[i,j] = attr[i]}.
#'
#' \code{"column"}(duplicating column actors' attributes: receiver effect): \eqn{Matrix[i,j] = attr[j]}.
#'
#' \code{"max"}(maximam value): \eqn{Matrix[i,j] = max(attr[i], attr[j])}.
#'
#' \code{"min"}(minimum value): \eqn{Matrix[i,j] = min(attr[i], attr[j])}.
#'
#' \code{"mean"}(mean value): \eqn{Matrix[i,j] = mean(attr[i], attr[j])}.
#'
#' @param attr vector of input attribute.
#' @param method \code{"match"}, \code{"match.one"}, \code{"diff"},
#'   \code{"abs.diff"}, \code{"sqd.diff"}, \code{"sum"}, \code{"product"},
#'   \code{"identity"}, \code{"row"}, \code{"column"}, \code{"max"},
#'   \code{"min"} or \code{"mean"}. Default is \code{"match"}.
#' @param diag.val the diag value to be filled. If \code{NULL} is set, no value
#'   is inputed.
#' @return a matrix converted from a vector.
#'
#' @examples
#' vec <- c(1,1,2,2,3,3)
#' names(vec) <- letters[1:6]
#' attr2mat(vec)
#' attr2mat(vec, "match.one")
#' attr2mat(vec, "match.one", diag.val = NULL)

attr2mat <- function(attr,
                     method = "match",
                     diag.val = 0){
  len <- length(attr)
  mat <- matrix(NA, len, len)

  # names
  dimnames(mat) <- list(names(attr), names(attr))

  # method
  if(method=="match"){
    for(i in 1:len){
      for(j in 1:len){
        if(attr[i]==attr[j]){
          mat[i,j] <- 1
        }else{
          mat[i,j] <- 0
        }
      }
    }
  }else if(method=="match.one"){
    for(i in 1:len){
      for(j in 1:len){
        if(attr[i]==1 & attr[j]==1){
          mat[i,j] <- 1
        }else{
          mat[i,j] <- 0
        }
      }
    }
  }else if(method=="diff"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- attr[i] - attr[j]
      }
    }
  }else if(method == "abs.diff"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- abs(attr[i] - attr[j])
      }
    }
  }else if(method=="sqd.diff"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- (attr[i]-attr[j])^2
      }
    }
  }else if(method=="sum"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- attr[i]+attr[j]
      }
    }
  }else if(method=="product"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- attr[i]*attr[j]
      }
    }
  }else if(method=="identity"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- 2*(attr[i]*attr[j])/(attr[i]^2+attr[j]^2)
      }
    }
  }else if(method=="row"){
    for(i in 1:len){
      mat[i,] <- attr[i]
    }
  }else if(method=="column"){
    for(j in 1:len){
      mat[,j] <- attr[j]
    }
  }else if(method=="max"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- max(c(attr[i],attr[j]))
      }
    }
  }else if(method=="min"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- min(c(attr[i], attr[j]))
      }
    }
  }else if(method=="mean"){
    for(i in 1:len){
      for(j in 1:len){
        mat[i,j] <- mean(c(attr[i],attr[j]))
      }
    }
  }else{
    stop("no method available")
  }

  # diag.val
  if(is.null(diag.val)){
    mat <- mat
  }else{
    diag(mat) <- diag.val
  }

  # return
  return(mat)
}
