#' Probabilities of maximization, by Beta PERT distributions
#' @description This function computes the probabilities of each alternative maximizing the preference per criterion, using Beta PERT distributions to randomize the decision matrix.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
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
#' s = 4 # Shape
#' PMax.Beta(x,s)
#' @importFrom  mc2d dpert ppert
#' @export
PMax.Beta = function (x,s) {

  m = x
  PMax = x

  max = apply(x,2,max)
  min = apply(x,2,min)

  for (j in 1:ncol(x))
  {
    for (i in 1:nrow(x))
    {
      PMax[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],m[,j][-i],max[j],s))*dpert(x,min[j],m[,j][[i]],max[j],s)}),min[j],max[j])) $value
    }}
  PMax[,]
  rownames(PMax) = paste0("Alt",1:nrow(m))
  colnames(PMax) = paste0("Crit",1:ncol(m))

  Result = PMax
  Result
}
