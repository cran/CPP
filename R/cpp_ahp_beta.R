#' CPP Additive Weighting with Probabilistic AHP using Beta PERT distributions
#' @description This function computes CPP by additive weighting. Experts' estimatives are based on pair-wise comparisons of criteria and are joined in a list of matrices. The estimatives are used as parameters to probabilistic distributions. The minimum, mean, and maximum values of each pair of criteria are used to model Beta PERT distributions. Randomic values are generated and applied to the AHP method. The matrix that comprises de minimum AHP Consistent Index is used to return the criteria weights.
#' @param n Random numbers created from Beta PERT distributions, using the parameters 'min', 'mean' and 'max' of each pair-wise criteria comparison elicited from the experts.
#' @param s Shape of a Beta PERT distribution, as described in Package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data elicited from experts.
#' @param list Pair-wise comparison matrices of expert opinions. The function 'list' is embedded in R.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @return Weights returned from the AHP method. PMax are the joint probabilities of each alternative being higher than the others, per criterion. CPP gives the final scores and ranks of alternatives by weighted sum.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @references Saaty, Thomas L. (1980). The analytic hierarchy process: planning, priority setting, resource allocation, McGraw-Hill.
#' @examples
#' n=5000 # simulation
#' s=6 # shape of Beta PERT distribution
#' # Expert pair-wise evaluations
#' Exp.1 = matrix(c(1,0.2,0.3,5,1,0.2,3,5,1),3,3)
#' Exp.2 = matrix(c(1,2,8,0.5,1,6,0.12,0.16,1),3,3)
#' Exp.3 = matrix(c(1,0.5,0.5,2,1,6,2,0.16,1),3,3)
#' Exp.4 = matrix(c(1,3,4,0.3,1,0.5,0.25,0.3,1),3,3)
#' Exp.5 = matrix(c(1,4,5,0.25,1,1,0.2,1,1),3,3)
#' list = list(Exp.1,Exp.2,Exp.3,Exp.4,Exp.5)
#' # Alternatives' original scores
#' Alt.1 = c(30,86,-5)
#' Alt.2 = c(26,77,-12)
#' Alt.3 = c(22,93,-4)
#' Alt.4 = c(34,65,-10)
#' Alt.5 = c(31,80,-8)
#' Alt.6 = c(29,79,-9)
#' Alt.7 = c(37,55,-15)
#' Alt.8 = c(21,69,-11)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8) # Decision matrix
#' CPP.AHP.Beta(n,s,list,x)
#' @importFrom mc2d dpert ppert
#' @importFrom stats integrate dnorm pnorm
#'@export
CPP.AHP.Beta = function (n,s,list,x){

  min = apply(simplify2array(list), 1:2, min)
  mean = apply(simplify2array(list), 1:2, mean)
  max = apply(simplify2array(list), 1:2, max)

  c=nrow(min)
  d = c^2
  simu = vector("list", d)
  k=1
  for (i in 1:c)
  {
    for (j in 1:c)
    {
      simu[[k]] <- rpert(n,min[i,j],mean[i,j],max[i,j],s)
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

  ### Saaty's Random Indices (RI)

  RI = c(0,0,0.58,0.9,1.12,1.24,1.32,1.41,1.45,1.49)

  min = which.min(CI)
  index = CI[[min]]/RI[c]
  w.min = weight[[min]]

  w = w.min

  ### Decision matrix normalization

  y = t(as.matrix(apply(x,2,sum)))
  dadosn=x
  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      dadosn[i,j] = x[i,j]/y[j]
    }}
  dadosn = replace(dadosn, dadosn == 0, 0.0000000001)

  ### PMax by Normal distributions

  x = dadosn
  PMax = x
  mat = x

  sd = apply(x,2,sd)

  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      PMax[i,j] = (integrate(Vectorize(function(x) {prod(pnorm(x,mat[,j][-i],sd[j]))*dnorm(x,mat[,j][[i]],sd[j])}),-2,2)) $value
    }}
  PMax = PMax[,]

  ### SAW

  saw = PMax%*%w
  rank = rank(-saw)

  SAW = cbind(saw, rank)
  colnames(SAW) = c("SAW","Rank")

  Result = list(Weights.AHP=w, PMax=PMax, CPP=SAW)
  Result
}
