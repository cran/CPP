#' CPP by axes using Beta PERT distributions
#' @description This function computes the CPP by axes, using Beta PERT distributions to randomize the decision matrix. The CPP by axes is used to rank alternatives in multicriteria decision problems. The "Progressive-Conservative" and the "Optimist-Pessimist" axes emulate four decision maker's points of view.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @return PMax are the joint probabilities of each alternative being higher than the others, per criterion. PMin are the joint probabilities of each alternative being lower than the others, also by criterion. Axes returns the alternatives' scores by axis and ranking for decisionmaking.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @references Garcia, Pauli A. A. & Sant'Anna, Annibal P. (2015). Vendor and logistics provider selection in the construction sector: A probabilistic preferences composition approach. Pesquisa Operacional 35.2: 363-375.
#' @examples
#' # Alternatives' original scores
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5) # Decision matrix
#' s = 4 # Shape
#' CPP.Axes.Beta(x,s)
#' @importFrom mc2d dpert ppert
#' @export
CPP.Axes.Beta = function (x,s) {

  PMax = x
  m = x

  max = apply(x,2,max)
  min = apply(x,2,min)

  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      PMax[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],m[,j][-i],max[j],shape=s))*dpert(x,min[j],m[,j][[i]],max[j],shape=s)}),min[j],max[j])) $value
    }}
  PMax = PMax[,]

  ####

  PMin = x

  max = apply(x,2,max)
  min = apply(x,2,min)

  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      PMin[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min[j],m[,j][-i],max[j],shape=s))*dpert(x,min[j],m[,j][[i]],max[j],shape=s)}),min[j],max[j])) $value
    }}
  PMin = PMin[,]

  ### Composition by axes

  # PP point of view

  PP = apply(PMax,1,prod)
  PP.rank = rank(-PP)

  # PO point of view

  Probs.m = 1-PMax
  PO = 1-(apply(Probs.m,1,prod))
  PO.rank = rank(-PO)

  # CP point of view

  Probs.mm = 1-PMin
  CP = apply(Probs.mm,1,prod)
  CP.rank = rank(-CP)

  # CO point of view

  CO = 1-(apply(PMin,1,prod))
  CO.rank = rank(-CO)

  Result = cbind(PP,PP.rank,PO,PO.rank,CP,CP.rank,CO,CO.rank)
  colnames(Result) = c("PP","R","PO","R","CP","R","CO","R")

  Result <- list(PMax=PMax, PMin=PMin, Axes=Result)
  Result
}
