#' CPP by Choquet integrals, using Beta PERT distributions
#' @description This function computes the CPP by Choquet integrals, using Beta PERT distributions to randomize the decision matrix. The CPP by Choquet integrals is used to rank alternatives in multicriteria decision problems.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data elicited from experts.
#' @return PMax are the joint probabilities of each alternative being higher than the others, per criterion. Capacities are the interactions of all combined criteria, computed by the Progressive-Optimistic (PO) point of view. Choq returns the alternatives' scores by Choquet integrals and their respetive rankings for decisionmaking. Shap returns the Shapley indices, which are associated with criteria weights.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @examples
#' # Alternatives' original scores
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5) # Decision matrix
#' s = 4 # Shape
#' CPP.Choquet.Beta(x,s)
#' @importFrom  kappalab capacity Choquet.integral Shapley.value
#' @importFrom  mc2d dpert ppert rpert
#' @importFrom  stats quantile
#' @importFrom  utils combn
#' @export
CPP.Choquet.Beta = function(x,s) {

  ### PMax Beta

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


  ### Capacities

  PMax.m = 1-PMax


  ### Union of subsets

  n = ncol(m)-1
  p = ncol(m)

  # Combinations 2:n

  syn = vector("list", n)
  for (i in 2:n)
  {
    syn[[i]] = 1-(t(apply(PMax.m, 1, function(x){combn(x,i,prod)})))
  }
  syn = do.call(cbind, syn)

  # Subset of all criteria
  ult.syn = 1-(apply(PMax.m, 1, function(x){combn(x,p,prod)}))

  ### MAX value determination

  mat.syn = cbind(PMax,syn,ult.syn)
  colnames(mat.syn)=NULL
  max1.syn = apply(mat.syn,2,max)
  max2.syn = max(max1.syn)

  ### capacities

  cap.syn = max1.syn/max2.syn
  kappa = c(0,cap.syn)
  kappa2= capacity(kappa)

  ### Choquet integrals and Shapley indices

  choquet = matrix(0,nrow(m),1)
  for (i in 1:nrow(m))
  {
    choquet[i] = Choquet.integral(kappa2,PMax[i,])
  }

  choquet.r = apply(-choquet,2,rank)
  choquet = cbind(choquet,choquet.r)
  rownames(choquet) = paste0("Alt",1:nrow(m))
  colnames(choquet) = c("CPP.Choquet","Rank")

  shapley = t(as.matrix(Shapley.value(kappa2)))
  colnames(shapley) = paste0("Crit",1:ncol(m))

  Result = list(PMax = PMax, Capacities = kappa2, Choq = choquet, Shap = shapley)
  Result
}
