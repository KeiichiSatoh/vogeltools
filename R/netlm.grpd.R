#' Linear Regression for Network Data with grouped random permutation
#'
#' \code{netlm} regresses the network variable in y on the network variables in
#' stack x using ordinary least squares. The resulting fits (and coefficients)
#' are then tested against the indicated null hypothesis.This is a extention of
#' \code{\link[sna]{netlm}}, as it supports the grouped network permutation.
#' This function is still in the testing phase.
#'
#' This is the extention of the \code{\link[sna]{netlm}} in the package
#' \code{sna}. This function supports the grouped network permutation when
#' calculating the significance level, although currently it support only the
#' method \code{qap} and \code{qapspp} (they are identical) for such
#' permutation. This function helps for QAP test in which the off-diagonal
#' matrix is set to be \code{NA}, for example.
#'
#' The tests in the \code{nullhyp} supported by \code{netlm.grpd} are
#' (currently) as follows:
#'
#' \code{classical}: tests based on classical asymptotics.
#'
#' \code{qap}: QAP permutation test; currently identical to \code{qapspp}.
#'
#' \code{qapspp}: QAP permutation test, using Dekker's "semi-partialling plus"
#' procedure.
#'
#' @param y dependent network variable. This should be a matrix, for obvious
#'   reasons; NAs are allowed, but dichotomous data is strongly discouraged due
#'   to the assumptions of the analysis.
#' @param x stack of independent network variables. Note that NAs are permitted,
#'   as is dichotomous data.
#' @param intercept logical; should an intercept term be added?
#' @param mode string indicating the type of graph being evaluated. "digraph"
#'   indicates that edges should be interpreted as directed; "graph" indicates
#'   that edges are undirected. mode is set to "digraph" by default.
#' @param diag logical; should the diagonal be treated as valid data? Set this
#'   true if and only if the data can contain loops. \code{diag} is \code{FALSE}
#'   by default.
#' @param grp.ns a vector indicates the number of actors in each group. The
#'   default is \code{NULL}. A raison d'etre of this extention function.
#' @param nullhyp string indicating the particular null hypothesis against which
#'   to test the observed estimands.
#' @param test.statistic string indicating the test statistic to be used for the
#'   Monte Carlo procedures.
#' @param tol tolerance parameter for \code{\link[base]{qr.solve}}.
#' @param reps integer indicating the number of draws to use for quantile
#'   estimation. (Relevant to the null hypothesis test only - the analysis
#'   itself is unaffected by this parameter.) Note that, as for all Monte Carlo
#'   procedures, convergence is slower for more extreme quantiles. By default,
#'   reps=1000.
#' @return An object of class \code{\link[sna]{netlm}}.
#'
#' @examples
#' # Example
#' library(sna)
#' xmat1.1 <- sna::rgraph(100)
#' xmat1.2 <- sna::rgraph(120)
#' xmat1 <- list2mat(list(xmat1.1,xmat1.2), fill = NA)
#'
#' xmat2.1 <- sna::rgraph(100)
#' xmat2.2 <- sna::rgraph(120)
#' xmat2 <- list2mat(list(xmat2.1,xmat2.2), fill = NA)
#'
#' mat.y1 <- sna::rgraph(100)
#' mat.y2 <- sna::rgraph(120)
#' ymat <- list2mat(list(mat.y1,mat.y2), fill = NA)
#'
#' res2.1 <- netlm.grpd(y = ymat, x = list(xmat1, xmat2),
#'                     grp.ns = c(100,120),reps = 100, nullhyp = "qap")
#' summary(res2.1)
#' res2.2 <- netlm.grpd(y = ymat, x = list(xmat1, xmat2),
#'                      grp.ns = c(100,120),nullhyp = "classical")
#' summary(res2.2)
#'
#' @export
#'
netlm.grpd <- function(
  y, x, intercept = TRUE,
  mode = "digraph", diag = FALSE,
  grp.ns = NULL,
  nullhyp = c("qap", "qapspp", "qapy", "qapx",
              "qapallx", "cugtie", "cugden", "cuguman", "classical"),
  test.statistic = c("t-value", "beta"), tol = 1e-07, reps = 1000){
  gettval <- function(x, y, tol) {
    xqr <- qr(x, tol = tol)
    coef <- qr.coef(xqr, y)
    resid <- qr.resid(xqr, y)
    rank <- xqr$rank
    n <- length(y)
    rdf <- n - rank
    resvar <- sum(resid^2)/rdf
    cvm <- chol2inv(xqr$qr)
    se <- sqrt(diag(cvm) * resvar)
    coef/se
  }
  gfit <- function(glist, mode, diag, tol, rety, tstat) {
    y <- sna::gvectorize(glist[[1]], mode = mode, diag = diag,
                    censor.as.na = TRUE)
    x <- vector()
    for (i in 2:length(glist)) x <- cbind(x, gvectorize(glist[[i]],
                                                        mode = mode, diag = diag, censor.as.na = TRUE))
    if (!is.matrix(x))
      x <- matrix(x, ncol = 1)
    mis <- is.na(y) | apply(is.na(x), 1, any)
    if (!rety) {
      if (tstat == "beta")
        qr.solve(x[!mis, ], y[!mis], tol = tol)
      else if (tstat == "t-value") {
        gettval(x[!mis, ], y[!mis], tol = tol)
      }
    }
    else {
      list(qr(x[!mis, ], tol = tol), y[!mis])
    }
  }
  y <- sna::as.sociomatrix.sna(y)
  x <- sna::as.sociomatrix.sna(x)
  if (is.list(y) || ((length(dim(y)) > 2) && (dim(y)[1] > 1)))
    stop("y must be a single graph in netlm.")
  if (length(dim(y)) > 2)
    y <- y[1, , ]
  if (is.list(x) || (dim(x)[2] != dim(y)[2]))
    stop("Homogeneous graph orders required in netlm.")
  nx <- stackcount(x) + intercept
  n <- dim(y)[2]
  g <- list(y)
  if (intercept)
    g[[2]] <- matrix(1, n, n)
  if (nx - intercept == 1)
    g[[2 + intercept]] <- x
  else for (i in 1:(nx - intercept)) g[[i + 1 + intercept]] <- x[i,
                                                                 , ]
  if (any(sapply(lapply(g, is.na), any)))
    warning("Missing data supplied to netlm; this may pose problems for certain null hypotheses.  Hope you know what you're doing....")
  fit.base <- gfit(g, mode = mode, diag = diag, tol = tol,
                   rety = TRUE)
  fit <- list()
  fit$coefficients <- qr.coef(fit.base[[1]], fit.base[[2]])
  fit$fitted.values <- qr.fitted(fit.base[[1]], fit.base[[2]])
  fit$residuals <- qr.resid(fit.base[[1]], fit.base[[2]])
  fit$qr <- fit.base[[1]]
  fit$rank <- fit.base[[1]]$rank
  fit$n <- length(fit.base[[2]])
  fit$df.residual <- fit$n - fit$rank
  tstat <- match.arg(test.statistic)
  if (tstat == "beta")
    fit$tstat <- fit$coefficients
  else if (tstat == "t-value")
    fit$tstat <- fit$coefficients/sqrt(diag(chol2inv(fit$qr$qr)) *
                                         sum(fit$residuals^2)/(fit$n - fit$rank))
  nullhyp <- match.arg(nullhyp)
  if ((nullhyp %in% c("qap", "qapspp")) && (nx ==
                                            1))
    nullhyp <- "qapy"
  if (nullhyp == "classical") {
    resvar <- sum(fit$residuals^2)/fit$df.residual
    cvm <- chol2inv(fit$qr$qr)
    se <- sqrt(diag(cvm) * resvar)
    tval <- fit$coefficients/se
    fit$dist <- NULL
    fit$pleeq <- pt(tval, fit$df.residual)
    fit$pgreq <- pt(tval, fit$df.residual, lower.tail = FALSE)
    fit$pgreqabs <- 2 * pt(abs(tval), fit$df.residual, lower.tail = FALSE)
  }
  else if (nullhyp %in% c("cugtie", "cugden", "cuguman")) {
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      gr <- g
      for (j in 1:reps) {
        gr[[i + 1]] <- switch(nullhyp, cugtie = rgraph(n,
                                                       mode = mode, diag = diag, replace = FALSE,
                                                       tielist = g[[i + 1]]), cugden = rgraph(n, tprob = gden(g[[i +
                                                                                                                   1]], mode = mode, diag = diag), mode = mode,
                                                                                              diag = diag), cuguman = (function(dc, n) {
                                                                                                rguman(1, n, mut = dc[1], asym = dc[2], null = dc[3],
                                                                                                       method = "exact")
                                                                                              })(dyad.census(g[[i + 1]]), n))
        repdist[j, i] <- gfit(gr, mode = mode, diag = diag,
                              tol = tol, rety = FALSE, tstat = tstat)[i]
      }
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="),
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="),
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat),
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapy") {
    repdist <- matrix(0, reps, nx)
    gr <- g
    for (i in 1:reps) {
      gr[[1]] <- rmperm.grpd(g[[1]], grp.ns)
      repdist[i, ] <- gfit(gr, mode = mode, diag = diag,
                           tol = tol, rety = FALSE, tstat = tstat)
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="),
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="),
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat),
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapx") {
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      gr <- g
      for (j in 1:reps) {
        gr[[i + 1]] <- rmperm.grpd(gr[[i + 1]], grp.ns)
        repdist[j, i] <- gfit(gr, mode = mode, diag = diag,
                              tol = tol, rety = FALSE, tstat = tstat)[i]
      }
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="),
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="),
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat),
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapallx") {
    repdist <- matrix(0, reps, nx)
    gr <- g
    for (i in 1:reps) {
      for (j in 1:nx) gr[[1 + j]] <- rmperm.grpd(g[[1 + j]], grp.ns)
      repdist[i, ] <- gfit(gr, mode = mode, diag = diag,
                           tol = tol, rety = FALSE, tstat = tstat)
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="),
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="),
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat),
                                ">="), 2, mean)
  }
  else if ((nullhyp == "qap") || (nullhyp == "qapspp")) {
    xsel <- matrix(TRUE, n, n)
    if (!diag){
      diag(xsel) <- FALSE}
    if (mode == "graph"){
      xsel[upper.tri(xsel)] <- FALSE}
    # put the information of NA in y matrix
    xsel[is.na(y)]<- FALSE
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      xfit <- gfit(g[1 + c(i, (1:nx)[-i])], mode = mode,
                   diag = diag, tol = tol, rety = TRUE, tstat = tstat)
      xres <- g[[1 + i]]
      xres[xsel] <- qr.resid(xfit[[1]], xfit[[2]])
      if (mode == "graph")
        xres[upper.tri(xres)] <- t(xres)[upper.tri(xres)]
      for (j in 1:reps) repdist[j, i] <- gfit(c(g[-(1 +
                                                      i)], list(rmperm.grpd(xres, grp.ns))), mode = mode, diag = diag,
                                              tol = tol, rety = FALSE, tstat = tstat)[nx]
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="),
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="),
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat),
                                ">="), 2, mean)
  }
  fit$nullhyp <- nullhyp
  fit$names <- paste("x", 1:(nx - intercept), sep = "")
  if (intercept)
    fit$names <- c("(intercept)", fit$names)
  fit$intercept <- intercept
  class(fit) <- "netlm"
  fit
}


