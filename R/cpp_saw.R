#' CPP by weighted sum, with weights informed by the user
#' @description This function computes the CPP-SAW, using Normal distributions and weights defined by the decision maker. The CPP-SAW is used to evaluate alternatives by weighted sum.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param w Vector of weights assigned by the decision maker. Weights are normalized, just in case their sum differs from the unity.
#' @return Weights repeat the parameter 'w' if sum the unity, otherwise are normalized. PMax indicates the joint probabilities of each alternative being higher than the others, per criterion. CPP returns the alternatives' scores by weighted sum, indicating the preference ranks for decisionmaking.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @examples
#' # Decision matrix
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' Alt.6 = c(6,29,79,-9)
#' Alt.7 = c(8,37,55,-15)
#' Alt.8 = c(10,21,69,-11)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8)
#' w = c(0.2,0.3,0.4,0.1)
#' CPP.SAW(x,w)
#' @export
CPP.SAW = function (x,w){

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
  apply(dadosn,2,sum)

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

  z = sum(w)
  w = w/z # Normalization of weights

  saw = PMax%*%w
  rank = rank(-saw)

  SAW = cbind(saw, rank)
  colnames(SAW) = c("CPP.SAW","Rank")

  Result = list(Weights=w, PMax=PMax, CPP=SAW)
  Result
}
