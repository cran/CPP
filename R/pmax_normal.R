#' Probabilities of maximization, by Normal distributions
#' @description This function computes the probabilities of each alternative maximizing the preference per criterion, using Normal distributions to randomize the decision matrix.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @return PMax are the joint probabilities of each alternative being higher than the others, per criterion.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @examples
#' # Decision matrix
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5)
#' PMax.Normal(x)
#' @export
PMax.Normal = function (x) {

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

  # PMax

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
  PMax[,]
  rownames(PMax) = paste0("Alt",1:nrow(mat))
  colnames(PMax) = paste0("Crit",1:ncol(mat))

  Result = PMax
  Result
}
