#' CPP for sorting alternatives in ordinal classes
#' @description This function computes the CPP-Tri, using Beta PERT distributions to randomize the decision matrix. The CPP-Tri is used to classify alternatives, indicating the order of classes, whose quantity is defined by the decision maker. The probabilities of each alternative being higher and lower than the classes' profiles are composed by the Progressive-Pessimist (PP) point of view.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param q Vector of quantiles, indicating the classes' profiles.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @return Prob.Plus are the probabilities of each alternative being higher than the classes'profiles. Prob.Minus are the probabilities of each alternative being lower than the classes'profiles. CPP.Tri returns the alternatives' classes.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @references Sant'Anna, Annibal P.; Costa, Helder G.; Pereira, Valdecy (2015). CPP-TRI: a sorting method based on the probabilistic composition of preferences. International Journal of Information and Decision Sciences 7.3, 193-212.
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
#' q = c(0.65,0.35) # quantiles of classes' profiles.
#' s = 4 # Shape
#' CPP.Tri.Beta(x,q,s)
#' @importFrom  mc2d dpert ppert
#' @export
CPP.Tri.Beta = function(x,q,s) {

  A = x
  min = apply(A,2,min)
  max = apply(A,2,max)

  Perfil = matrix(0,length(q),ncol(A))

  for (a in 1:length(q))
  {
    Perfil[a,] = apply(A,2,function(x) quantile(x,q[a]))
  }

  A = rbind(Perfil,A)

  ##### Prob.Max OVER profiles

  listac = vector("list", length(q))

  for (a in 1:length(q))
  {
    listac[[a]] = matrix(0,nrow(A),ncol(A))

    for (j in 1:ncol(A))
    {
      for (i in 1:nrow (A))
      {
        listac[[a]][i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],A[a,j],max[j],s))*dpert(x,min[j],A[,j][[i]],max[j],s)}),min[j],max[j]))$value
      }}}
  names(listac) = paste0("Profile",1:length(q))
  listac

  ##### Prob.Min UNDER profiles

  listab = vector("list", length(q))

  for (a in 1:length(q))
  {
    listab[[a]] = matrix(0,nrow(A),ncol(A))

    for (j in 1:ncol(A))
    {
      for (i in 1:nrow (A))
      {
        listab[[a]][i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min[j],A[a,j],max[j],s))*dpert(x,min[j],A[,j][[i]],max[j],s)}),min[j],max[j]))$value
      }}}
  names(listab) = paste0("Profile",1:length(q))
  listab

  ###################
  ##### Sorting #####
  ###################

  PPacima = lapply(listac,apply,1,prod)
  PPabaixo = lapply(listab,apply,1,prod)

  Difs = vector("list", length(q))
  for (a in 1:length(q))
  {
    Difs[[a]] = abs(PPacima[[a]]-PPabaixo[[a]])
  }

  Difs = t(do.call(rbind, Difs))

  Classe = as.matrix(apply(Difs,1,which.min))
  Classe = Classe[-c(1:length(q)),]
  Classe = matrix(Classe,nrow(x),1)
  rownames(Classe) = paste0("Alt",1:nrow(x))
  colnames(Classe) = c("Class")

  ### Results

  PPacima = t(do.call(rbind, PPacima))
  PPacima = PPacima[-c(1:length(q)),]

  PPabaixo = t(do.call(rbind, PPabaixo))
  PPabaixo = PPabaixo[-c(1:length(q)),]

  Result = list(Prob.Plus = PPacima,Prob.Minus = PPabaixo, CPP.Tri = Classe)
  Result
}
