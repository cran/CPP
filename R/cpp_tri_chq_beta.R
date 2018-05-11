#' CPP for sorting alternatives, based on Choquet integrals
#' @description This function computes the CPP-Tri with Choquet integrals, using Beta PERT distributions to randomize the decision matrix. The CPP Tri is used to classify alternatives, indicating the order of a number of classes defined by the decision maker. The probabilities of each alternative being higher and lower than the classes' profiles are composed by Choquet integrals.
#' @param x Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria negative.
#' @param q Vector of quantiles, indicating the classes' profiles.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @return Choquet.O returns the probabilities of each alternative being higher than the classes' profiles, composed by Choquet integrals. Choquet.U indicates the probabilities of each alternative being lower than the classes' profiles. CPP.Tri.Chq returns the alternatives' classes.
#' @references Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
#' @references Sant’Anna, Annibal P.; Lima, Gilson B. A.; Gavião, Luiz O. (2018) A Probabilistic approach to the inequality adjustment of the Human Development Index. Pesquisa Operacional 38.1: 99-116.
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
#' CPP.Tri.Choquet(x,q,s)
#' @importFrom  kappalab capacity Choquet.integral Shapley.value
#' @importFrom  mc2d dpert ppert
#' @export
CPP.Tri.Choquet = function(x,q,s){

  ############
  ### PMax ###
  ############

  m = x
  PMax = m

  max = apply(m,2,max)
  min = apply(m,2,min)

  for (j in 1:ncol(m))
  {
    for (i in 1:nrow(m))
    {
      PMax[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],m[,j][-i],max[j],s))*dpert(x,min[j],m[,j][[i]],max[j],s)}),min[j],max[j])) $value
    }}

  ##############################
  ### Capacity determination ###
  ##############################

  PMax.m = 1-PMax

  n = ncol(m)-1
  p = ncol(m)

  syn = vector("list", n)
  for (i in 2:n)
  {
    syn[[i]] = 1-(t(apply(PMax.m, 1, function(x){combn(x,i,prod)})))
  }
  syn = do.call(cbind, syn)

  # Last criteria
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


  #########################
  ### Sorting procedure ###
  #########################

  A = x
  min = apply(A,2,min)
  max = apply(A,2,max)

  Perfil = matrix(0,length(q),ncol(A))

  for (a in 1:length(q))
  {
    Perfil[a,] = apply(m,2,function(x) quantile(x,q[a]))
  }

  A = rbind(Perfil,A)


  ##### Probs OVER profiles

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


  ##### Probs UNDER profiles

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


  #########################
  ### Choquet integrals ###
  #########################

  choquet.ac = vector("list", length(q))

  for (b in 1:length(q))
  {
    for (i in 1:nrow(x))
    {
      choquet.ac[[b]][i] = Choquet.integral(kappa2,listac[[b]][i,])
    }}

  choquet.ab = vector("list", length(q))

  for (b in 1:length(q))
  {
    for (i in 1:nrow(x))
    {
      choquet.ab[[b]][i] = Choquet.integral(kappa2,listab[[b]][i,])
    }}


  Difs = vector("list", length(q))
  for (b in 1:length(q))
  {
    Difs[[b]] = abs(choquet.ac[[b]]-choquet.ab[[b]])
  }

  Difs = t(do.call(rbind, Difs))

  Classe = as.matrix(apply(Difs,1,which.min))
  rownames(Classe) = paste0("Alt",1:nrow(x))
  colnames(Classe) = c("Class")

  Result = list(Choquet.O = choquet.ac, Choquet.U = choquet.ab, CPP.Tri.Chq = Classe)
  Result
}
