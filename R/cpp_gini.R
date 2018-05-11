#' CPP by the Gini Index, using Beta PERT distributions
#' @description The CPP by the Gini Index is used to rank alternatives by evenness of evaluations, in multicriteria decision problems.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis.
#' @return PMax are the joint probabilities of each alternative being higher than the others, per criterion. CPP.Gini returns the alternatives' scores by the Gini Index and their respective preference ranks for decisionmaking.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer
#' @references Gaviao, Luiz O. & Lima, Gilson B.A. (2017) Support decision to player selection: an application of the CPP in soccer, Novas Edições Acadêmicas [in Portuguese].
#' @examples
#' # Alternatives' original scores
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5) # Decision matrix
#' s = 4 # Shape
#' CPP.Gini(x,s)
#' @importFrom  ineq ineq
#' @importFrom  mc2d dpert ppert
#' @export
CPP.Gini = function(x,s) {

  b = x
  PMax = x

  max = apply(b,2,max)
  min = apply(b,2,min)

  for (j in 1:ncol(b))
  {
    for (i in 1:nrow(b))
    {
      PMax[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],b[,j][-i],max[j],s))*dpert(x,min[j],b[,j][[i]],max[j],s)}),min[j],max[j])) $value
    }}

  b = PMax[,]

  ### Gini Index

  gin = matrix(nrow(b),1)

  for (i in 1:nrow(b))
  {
    m = as.vector(b[i,])
    gin[i] = ineq(m, parameter = NULL, type = "Gini")
  }

  r.gin = rank(gin)
  index = cbind(gin,r.gin)
  colnames(index) = c("Index","Rank")
  Result = list(PMax = PMax,CPP.Gini = index)
  Result
}
