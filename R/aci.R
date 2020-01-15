#' Calculating the Advocacy Coalition Index.
#'
#' \code{aci} calculates the advocacy coalition index.
#'
#' This function returns the adovcacy coalition index (ACI) and homophily matrix.
#'
#'
#' @param matrix a matrix.
#' @param policy.score a vector indicating the policy score ranging 0 to 1.
#' @param cut.value the cut-off value to differentiate the homopilious ties. The
#'   default value is 0.5.
#' @return a list object containing,
#'
#' (1) actor level broker-ties, missing-ties;score as well as ACI (in-coming, and out-going score differentiated);
#'
#' (2) the whole network level ACI, broker ties score and missing-ties score;
#'
#' (3) coalition matrix;
#'
#' (4) original matrix supplied as matrix;
#'
#' (5) preference agreement matrix.
#' @references Satoh, Keiichi, Antti Gronow, and Tuomas Yla-Anttila. (2019). "The Advocacy Coalition Index: Identifying Coalitions by Simultaneously Taking into Account Coordination of Action and Belief Homophily". Paper Presented at ECPR 2019.(\url{https://ecpr.eu/Events/PaperDetails.aspx?PaperID=46846&EventID=123}).
#' @examples
#' # matrix
#' mat <- matrix(c(0,1,1,0,0,0,
#'                 1,0,0,0,0,0,
#'                 1,0,0,0,0,1,
#'                 0,0,0,0,1,1,
#'                 0,0,0,1,0,1,
#'                 0,0,1,1,1,0),6,6,byrow = T)
#'
#' namelist <- c(paste(rep("A",3),1:3, sep = ""),
#'               paste(rep("B",3),1:3, sep = ""))
#' dimnames(mat)<-list(namelist,namelist)
#'
#' # policy score
#' policy.score <- c(1,1,1,0,0,0)
#' names(policy.score) <- namelist
#'
#' # ACI
#' ACI <- aci(matrix = mat,                # imput matrix
#'              policy.score = policy.score, # imput policy score
#'              cut.value = 0.5)             # cut-off value (alpha) for the coalition matrix
#'
#' ACI$act                 # actor-level ACI
#' ACI$whole               # network-level  ACI
#' ACI$coalition.mat       # coalition matrix
#' ACI$original.matrix.Y
#' ACI$preference.agreement.matrix.X
#'
# Example 2: completely in line with ACF
#' mat2 <- mat
#' mat2[2,3]<-1;mat2[3,2]<-1
#' mat2[3,6]<-0;mat2[6,3]<-0
#'
#' aci(mat2,policy.score)
#'
#' Example 3: completely opposite to ACF
#' mat3 <- mat2
#' mat3 <- 1-mat3
#' aci(mat3,policy.score)
#'
#' @export
aci <- function(matrix,
                policy.score,
                cut.value = 0.5){
  # imput data
  Y <- matrix
  a <- policy.score
  alpha <- cut.value
  diag(Y)<- 0

  # number of actors
  n <- length(a)

  # transform the score into matrix
  X <- matrix(NA, n, n)
  dimnames(X) <- dimnames(Y)
  for(i in 1:length(a)){for(j in 1:length(a)){
    X[i,j] <- 1-abs(a[i]-a[j])
  }}
  diag(X) <- NA

  # distance(D)
  D <- abs(Y-X)


  # differentiate the submatrix
  D.broker <- D
  D.broker[Y<X] <- 0
  D.missing <- D
  D.missing[Y>X] <- 0

  # coalition matrix
  coalition.mat <- Y
  coalition.mat[D>alpha] <- 0
  dimnames(coalition.mat)<-dimnames(Y)

  # ACI-actor(out)
  act.broker.out <- rowSums(D.broker, na.rm = T)/rowSums(Y)
  act.missing.out  <- rowSums(D.missing, na.rm = T)/rowSums(X, na.rm = T)
  act.ACI.out <- 1-(rowSums(D, na.rm = T)/(n-1))
  # ACI-actor(in)
  act.broker.in <- colSums(D.broker, na.rm = T)/colSums(Y)
  act.missing.in  <- colSums(D.missing, na.rm = T)/colSums(X, na.rm = T)
  act.ACI.in    <- 1-(colSums(D, na.rm = T)/(n-1))
  # ACI-whole
  whole.broker  <- sum(D.broker, na.rm = T)/sum(Y)
  whole.missing   <- sum(D.missing, na.rm = T)/sum(X, na.rm = T)
  whole.ACI     <- 1-(sum(D, na.rm = T)/(n*(n-1)))
  ## attach th result
  act <- data.frame(act.ACI.out    = act.ACI.out,
                    act.broker.out = act.broker.out,
                    act.missing.out  = act.missing.out,
                    act.ACI.in     = act.ACI.in,
                    act.broker.in  = act.broker.in,
                    act.missing.in   = act.missing.in)
  whole <- c(whole.ACI     = whole.ACI,
             whole.broker  = whole.broker,
             whole.missing   = whole.missing)
  result <- list(act=act,
                 whole=whole,
                 coalition.mat = coalition.mat,
                 original.matrix.Y = Y,
                 preference.agreement.matrix.X = X)
  ## return the result
  print(result$whole)
  invisible(result)
}
