#' CPP with multiple perspectives for human resources evaluation
#' @description This function computes the CPP-rh, using Beta PERT distributions to randomize the decision matrices. The CPP-rh is used to evaluate alternatives by integrating the CPP-Tri, the CPP-Malmquist, the CPP-Gini and the CPP by axes. The CPP-rh and the CPP-mb are very similar, but the CPP-rh does not include the alternatives's market value.
#' @param t1 Decision matrix of Alternatives (rows) and Criteria (columns) in the previous moment '1'. Benefit criteria must be positive and cost criteria must be negative.
#' @param t2 Decision matrix of Alternatives (rows) and Criteria (columns) in the following moment '2'. Benefit criteria must be positive and cost criteria must be negative.
#' @param q Vector of quantiles, indicating the classes' profiles.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @return Class identifies the alternatives' classes, according to the selected profiles. CPP-RH returns the alternatives' scores per class.
#' @examples
#' ## Decision matrix of the previous moment '1'.
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' Alt.6 = c(6,29,79,-9)
#' Alt.7 = c(8,37,55,-15)
#' Alt.8 = c(10,21,69,-11)
#' t1 = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8)
#' ## Decision matrix of the following moment '2'.
#' Alt.1 = c(3,29,82,-3)
#' Alt.2 = c(6,28,70,-8)
#' Alt.3 = c(2,20,99,-8)
#' Alt.4 = c(5,31,62,-14)
#' Alt.5 = c(9,27,73,-5)
#' Alt.6 = c(4,33,85,-13)
#' Alt.7 = c(9,39,59,-10)
#' Alt.8 = c(8,19,77,-9)
#' t2 = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5,Alt.6,Alt.7,Alt.8)
#' q = c(0.65,0.35) # quantiles of class profiles
#' s = 4 # Shape
#' CPP.rh(t1,t2,q,s)
#' @importFrom  ineq ineq
#' @importFrom  mc2d dpert ppert
#' @export
CPP.rh = function(t1,t2,q,s){

  rownames(t1) = paste('Alt', 1:(nrow(t1)))
  rownames(t2) = paste('Alt', 1:(nrow(t1)))

  ### CPP Tri

  A = t2
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


  ### Sorting procedure

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

  ### Original matrices "t1" and "t2"

  t1 = cbind(t1,Classe)
  t2 = cbind(t2,Classe)
  c = ncol(t1)

  list1 = vector("list", length(q))
  list2 = vector("list", length(q))
  for (a in 1:length(q))
  {
    list1[[a]] = subset(t1[,-c(c)], t1[,c] == a)
    list2[[a]] = subset(t2[,-c(c)], t2[,c] == a)
  }


  ### CPP Gini

  PMax = vector("list", length(q))

  for (a in 1:length(q))
  {
    b = list2[[a]]
    PMax[[a]] = matrix(0,nrow(b),ncol(b))
    max = apply(b,2,max)
    min = apply(b,2,min)
    for (j in 1:ncol(b))
    {
      for (i in 1:nrow(b))
      {
        PMax[[a]][i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],b[,j][-i],max[j],s))*dpert(x,min[j],b[,j][[i]],max[j],s)}),min[j],max[j])) $value
      }}}


  gin = vector("list", length(q))
  for (a in 1:length(q))
  {
    b = PMax[[a]]
    for (i in 1:nrow(b))
    {
      g = b[i,]
      gin[[a]][i] = ineq(g, parameter = NULL, type = "Gini")
    }
    gin[[a]] = as.matrix(gin[[a]])
    rownames(gin[[a]]) = c(rownames(list1[[a]]))
  }


  ### CPP Malmquist

  c = ncol(t1)
  pre = t1[,-c(c)]
  pos = t2[,-c(c)]

  min.pre = apply(pre, 2, min) # Minimum value per criterion
  max.pre = apply(pre, 2, max) # Maximum value per criterion

  min.pos = apply(pos, 2, min) # Minimum value per criterion
  max.pos = apply(pos, 2, max) # Maximum value per criterion


  # Malmquist - Conservative (MC)

  MC_pre_t1 = pre
  MC_pre_t = pre
  MC_pos_t1 = pos
  MC_pos_t = pos

  # MC PRE

  # MC PRE T+1

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MC_pre_t1[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MC_pre_t1 = 1-(MC_pre_t1)
  MC_pre_t1 = apply(MC_pre_t1,1,prod)


  # MC PRE T

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MC_pre_t[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MC_pre_t = 1-(MC_pre_t)
  MC_pre_t = apply(MC_pre_t,1,prod)

  MC_PRE = (MC_pre_t1)/(MC_pre_t)


  # MC POS

  # MC POS T+1

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MC_pos_t1[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MC_pos_t1 = 1-(MC_pos_t1)
  MC_pos_t1 = apply(MC_pos_t1,1,prod)


  # MC POS T

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MC_pos_t[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MC_pos_t = 1-(MC_pos_t)
  MC_pos_t = apply(MC_pos_t,1,prod)


  MC_POS = (MC_pos_t1)/(MC_pos_t)


  # MC Index
  MC = sqrt((MC_PRE)*(MC_POS))


  # Malmquist - Progressive (MP)

  # Variables

  MP_pre_t1 = pre
  MP_pre_t = pre
  MP_pos_t1 = pos
  MP_pos_t = pos


  # MP PRE

  # MP PRE T+1

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MP_pre_t1[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MP_pre_t1 = 1-(MP_pre_t1)
  MP_pre_t1 = apply(MP_pre_t1,1,prod)


  # MP PRE T

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MP_pre_t[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],)}),min.pre[j],max.pre[j]))$value
    }}

  MP_pre_t = 1-(MP_pre_t)
  MP_pre_t = apply(MP_pre_t,1,prod)


  # Resultado MP PRE
  MP_PRE = (MP_pre_t1)/(MP_pre_t)


  # MP POS

  # MP POS T+1

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MP_pos_t1[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MP_pos_t1 = 1-(MP_pos_t1)
  MP_pos_t1 = apply(MP_pos_t1,1,prod)


  # MP POS T

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MP_pos_t[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MP_pos_t = 1-(MP_pos_t)
  MP_pos_t = apply(MP_pos_t,1,prod)

  MP_POS = (MP_pos_t1)/(MP_pos_t)


  # MP Index
  MP = sqrt((MP_PRE)*(MP_POS))


  # PROB. MALMQUIST INDEX (IMP)

  IMP = sqrt(MC/MP)

  IMP1 = cbind(IMP,Classe)
  c = ncol(IMP1)

  IMP2 = vector("list", length(q))
  for (a in 1:length(q))
  {
    IMP2[[a]] = subset(IMP1[,-c(c)], IMP1[,c] == a)
    #rownames(IMP2[[a]]) = c(rownames(list1[[a]]))
  }

  ### CPP RH

  matrix = vector("list", length(q))
  for (a in 1:length(q))
  {
    matrix[[a]] = cbind(-gin[[a]],IMP2[[a]])
  }

  for(i in seq_along(matrix))
  {
    colnames(matrix[[i]]) <- c("CPP.Gini","CPP.Malmquist")
  }

  PMax.rh = vector("list", length(q))
  for (a in 1:length(q))
  {
    b = matrix[[a]]
    PMax.rh[[a]] = matrix(0,nrow(b),ncol(b))
    max = c(0,max(b[,2]))
    min = c(-1,min(b[,2]))
    for (j in 1:ncol(b))
    {
      for (i in 1:nrow(b))
      {
        PMax.rh[[a]][i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],b[,j][-i],max[j],s))*dpert(x,min[j],b[,j][[i]],max[j],s)}),min[j],max[j])) $value
      }}}

  for(i in seq_along(PMax.rh))
  {
    colnames(PMax.rh[[i]]) <- c("CPP.Gini","CPP.Malmquist")
  }

  PP = vector("list", length(q))
  PP.r = vector("list", length(q))
  rh = vector("list", length(q))

  for (a in 1:length(q))
  {
    PP[[a]] = apply(PMax.rh[[a]],1,prod)
    PP.r[[a]] = rank(-PP[[a]])
    rh[[a]] = cbind(PMax.rh[[a]],PP[[a]],PP.r[[a]])
    rownames(rh[[a]]) = c(rownames(list1[[a]]))
  }

  for(i in seq_along(rh))
  {
    colnames(rh[[i]]) <- c("PMax.Gini","PMax.Malmquist","CPP.PP","Rank")
  }

  ### Results

  Results = list(Class = matrix, CPP.RH = rh)
  Results
}
