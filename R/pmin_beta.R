#' Probabilities of minimization, by Beta PERT distributions
#' @description This function computes the Probabilities of each alternative minimizing the preference per criterion, using Beta PERT distributions to randomize the decision matrix.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @return PMin are the joint probabilities of each alternative being lower than the others, per criterion.
#' @examples
#' # Decision matrix
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5)
#' s = 4 # Shape
#' PMin.Beta(x,s)
#' @importFrom  mc2d dpert ppert
#' @export
PMin.Beta = function (x,s) {

  m = x
  PMin = x

  max = apply(x,2,max)
  min = apply(x,2,min)

  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      PMin[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min[j],m[,j][-i],max[j],s))*dpert(x,min[j],m[,j][[i]],max[j],s)}),min[j],max[j])) $value
    }}
  PMin[,]
  rownames(PMin) = paste0("Alt",1:nrow(m))
  colnames(PMin) = paste0("Crit",1:ncol(m))

  Result = PMin
  Result
}
