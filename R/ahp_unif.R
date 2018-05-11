#' Probabilistic AHP using Uniform distributions
#' @description This function computes criteria weights, using AHP and randomic pair-wise evaluations by Uniform distributions.
#' @param n Random numbers created from Uniform distributions, using the parameters 'min' and 'max' of each pair-wise criteria comparison elicited from the experts.
#' @param list Pair-wise comparison matrices of expert opinions. The function 'list' is embedded in R.
#' @return Weights returned from a simulation of AHP with Uniform distributions. The weights are driven from the simulated matrix that gives the minimum AHP Consistent Index.
#' @references Saaty, Thomas L. (1980). The analytic hierarchy process: planning, priority setting, resource allocation, McGraw-Hill.
#' @examples
#' n=5000 # Simulation
#'# Expert pair-wise evaluations
#' Exp.1 = matrix(c(1,0.2,0.3,5,1,0.2,3,5,1),3,3)
#' Exp.2 = matrix(c(1,2,8,0.5,1,6,0.12,0.16,1),3,3)
#' Exp.3 = matrix(c(1,0.5,0.5,2,1,6,2,0.16,1),3,3)
#' Exp.4 = matrix(c(1,3,4,0.3,1,0.5,0.25,0.3,1),3,3)
#' Exp.5 = matrix(c(1,4,5,0.25,1,1,0.2,1,1),3,3)
#' list = list(Exp.1,Exp.2,Exp.3,Exp.4,Exp.5)
#' AHP.Unif(n,list)
#' @importFrom stats runif
#' @export
AHP.Unif = function (n,list) {

  min = apply(simplify2array(list), 1:2, min)
  max = apply(simplify2array(list), 1:2, max)

  c=nrow(min)
  d = c^2
  simu = vector("list", d)
  k=1
  for (i in 1:c)
  {
    for (j in 1:c)
    {
      simu[[k]] <- runif(n,min[i,j],max[i,j])
      k = k+1
    }}

  un = unlist(simu)
  abc = matrix(un,c^2,n,byrow=TRUE)

  ### AHP

  m = vector("list", n)
  weight = vector("list", n)
  CI = vector("list", n)

  for (a in 1:n)
  {
    m[[a]] = as.vector(abc[,a])
    matrix = matrix(m[[a]],nrow = c,ncol = c, byrow = TRUE)
    for (i in 1:c)
    {
      weight[[a]][i] <- prod(matrix[i,])^(1/c)
    }
    temp_sum <- sum(weight[[a]])
    weight[[a]] <- weight[[a]]/temp_sum
    lambda_max <- Re(eigen(matrix)$values[1])
    CI[[a]] <- (lambda_max-c)/(c-1)
  }

  ### Saaty's Random Index

  RI = c(0,0,0.58,0.9,1.12,1.24,1.32,1.41,1.45,1.49)

  min = which.min(CI)
  index = CI[[min]]/RI[c]
  w.min = weight[[min]]

  Result = c(Weight=w.min)
  Result
}
