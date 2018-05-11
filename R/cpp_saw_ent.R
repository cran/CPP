#' CPP by weighted sum, with weights computed from Shannon entropy.
#' @description This function computes the CPP-SAW, using Normal distributions to randomize the decision matrix and weights defined by entropy. The CPP-SAW Entropy is used to evaluate alternatives by weighted sum.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @return Weights by entropy.PMax are the joint probabilities of each alternative being higher than the others, per criterion. CPP returns the alternatives' scores by weighted sum, indicating the preference ranks for decisionmaking.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @examples
#' ## Decision matrix.
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' Alt.6 = c(6,29,79,-9)
#' Alt.7 = c(8,37,55,-15)
#' Alt.8 = c(10,21,69,-11)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8)
#' CPP.SAW.Entropy(x)
#' @export
CPP.SAW.Entropy = function (x){

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

  # Entropy per criteria (less dispersion indicates a higher entropy and a lower criteria weight)
  k = 1/log(nrow(x)) # constant
  fun = apply(dadosn, 2, function(x) x*log(x))
  ent = apply(fun, 2, sum)
  entropia = -k*ent

  # Diferential factor
  d = 1 - entropia

  # Normalization by sum
  w = as.matrix(d/(sum(d)))

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
  colnames(SAW) = c("CPP.SAW.Ent","Rank")

  Result = list(Weights=w, PMax=PMax, CPP=SAW)
  Result
}
