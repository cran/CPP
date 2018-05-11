#' Weights by entropy
#' @description This function computes weights by Shannon's entropy.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria negative.
#' @return Weights for each criterion.
#' @references Pomerol, Jean-Charles & Barba-Romero, Sergio. (2012) Multicriterion Decision in Management: Principles and Practice, Springer.
#' @examples
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' Alt.6 = c(6,29,79,-9)
#' Alt.7 = c(8,37,55,-15)
#' Alt.8 = c(10,21,69,-11)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8) # Decision matrix.
#' Entrop.weights(x)
#' @export
Entrop.weights = function (x) {

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
  w = d/(sum(d))
  w = c(Weights=w)
  w
}
