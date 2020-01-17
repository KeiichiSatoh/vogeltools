#' Calculate vif for network regression.
#'
#' \code{netvif} vectorizes the object supplied for the QAP regression and
#' passes to \code{vif} in the R package \code{car} to calculate the vif.
#'
#' This is a utility function for selecting the independent variables when
#' conducting QAP regression.
#'
#' @param y dependent network variable.
#' @param x stack of independent network variables (array or list object).
#' @param mode string indicating the type of graph being evaluated.
#'   \code{"digraph"}(default: directed graph) or \code{"graph"}(undirected
#'   graph) .
#' @param diag logical. Whether the diagonal value is considered or not
#'   (default).
#' @param fun intiger indicating the regression model: \code{"netlm"}(default)
#'   or \code{"netlogit"}.
#' @return a vector indicating vifs.
#' @seealso \code{\link[car]{vif}}, \code{\link[sna]{netlm}},
#'   \code{\link[sna]{netlogit}}.
#'
#' @examples
#' # Create some input graphs
#' x <- sna::rgraph(20,4)
#'
#' #Create a response structure
#' y <- x[1,,]+4*x[2,,]+2*x[3,,]   #Note that the fourth graph is unrelated
#'
#' @export

netvif <- function(y, x, mode = "digraph", diag = FALSE, fun = "netlm"){
  # vectorize the x-matrix
  x.vec <- sna::gvectorize(x, mode = mode, diag = diag)
  y.vec <- sna::gvectorize(y, mode = mode, diag = diag)
  # attach the names
  if(class(x)=="array"){
    colnames(x.vec) <- dimnames(x)[[1]]
  }else if(class(x)=="list"){
    colnames(x.vec) <- names(x)
  }else if(class(x)=="matrix"){
    colnames(x.vec) <- NULL
  }
  # attach the names is x.vec col is NULL
  if(is.null(colnames(x.vec))){
    colnames(x.vec) <- paste(rep("X",ncol(x.vec)),1:ncol(x.vec), sep = "")
  }
  # attach as data.frame
  dat <- data.frame(y.vec=y.vec,
                    data.frame(x.vec))

  # calculate the regression
  if(fun=="netlm"){
    mod <- lm(y.vec~., data = dat)
  }else if(fun=="netlogit"){
    mod <- glm(y.vec~., data = dat, family = binomial)
  }
  # pass to car::vif
  res <- car::vif(mod)
  # return
  return(res)
}
