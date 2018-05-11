#' CPP by the Malmquist Index, using Beta PERT distributions
#' @description The CPP-Malmquist is used to dynamic evaluation of alternatives, in multicriteria problems, considering two different moments.
#' @param m1 Decision matrix of Alternatives (rows) and Criteria (columns) in moment '1'. Benefit criteria must be positive and cost criteria must be negative.
#' @param m2 Decision matrix of Alternatives (rows) and Criteria (columns) in the following moment '2'. Benefit criteria must be positive and cost criteria must be negative.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis, which emulates the precision of data.
#' @return MC gives the Malmquist Conservative index. MP gives the Malmquist Progressive index. Finally, Index gives the CPP-Malmquist of all alternatives and their rankings for decisionmaking. The indices greater than one represent a relative evolution of the alternative between the two periods, while the indices lower than one reveal the alternatives that decreased performance in relation to the others.
#' @examples
#' # Alternatives' original scores
#' Alt.1 = c(2,30,86,-5)
#' Alt.2 = c(4,26,77,-12)
#' Alt.3 = c(3,22,93,-4)
#' Alt.4 = c(6,34,65,-10)
#' Alt.5 = c(5,31,80,-8)
#' m1 = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5) # Decision matrix of the previous moment '1'.
#' Alt.1 = c(3,29,82,-3)
#' Alt.2 = c(6,28,70,-8)
#' Alt.3 = c(2,20,99,-8)
#' Alt.4 = c(5,31,62,-14)
#' Alt.5 = c(9,27,73,-5)
#' m2 = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5) # Decision matrix of the following moment '2'.
#' s = 4 # Shape
#' CPP.Malmquist.Beta(m1,m2,s)
#' @importFrom  mc2d dpert ppert
#' @export
CPP.Malmquist.Beta = function(m1,m2,s) {

  pre = m1
  pos = m2

  ### Assumption of imprecision/uncertainty with Beta PERT Distributions
  ### Parameters: dpert(x,"min","moda","max","shape")


  ### MAX MIN values of moments "pre" and "pos"

  min.pre = apply(pre, 2, min) # Minimum value per criterion
  max.pre = apply(pre, 2, max) # Maximum value per criterion

  min.pos = apply(pos, 2, min) # Minimum value per criterion
  max.pos = apply(pos, 2, max) # Maximum value per criterion


  #################################
  # Malmquist - Conservative (MC) #
  #################################

  ### Variables

  MC_pre_t1 = pre
  MC_pre_t = pre
  MC_pos_t1 = pos
  MC_pos_t = pos

  ##############
  ### MC PRE ###
  ##############

  ### MC PRE T+1

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MC_pre_t1[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MC_pre_t1 = 1-(MC_pre_t1)
  MC_pre_t1 = apply(MC_pre_t1,1,prod)


  ### MC PRE T

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MC_pre_t[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MC_pre_t = 1-(MC_pre_t)
  MC_pre_t = apply(MC_pre_t,1,prod)


  ### Resultado MC PRE
  MC_PRE = (MC_pre_t1)/(MC_pre_t)


  ##############
  ### MC POS ###
  ##############

  ### MC POS T+1

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MC_pos_t1[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MC_pos_t1 = 1-(MC_pos_t1)
  MC_pos_t1 = apply(MC_pos_t1,1,prod)


  ### MC POS T

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MC_pos_t[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MC_pos_t = 1-(MC_pos_t)
  MC_pos_t = apply(MC_pos_t,1,prod)


  ### Resultado MC POS
  MC_POS = (MC_pos_t1)/(MC_pos_t)


  ##################
  #### MC Index ####
  MC = sqrt((MC_PRE)*(MC_POS))


  ################################
  # Malmquist - Progressive (MP) #
  ################################

  ### Variables

  MP_pre_t1 = pre
  MP_pre_t = pre
  MP_pos_t1 = pos
  MP_pos_t = pos

  ##############
  ### MP PRE ###
  ##############

  ### MP PRE T+1

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MP_pre_t1[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MP_pre_t1 = 1-(MP_pre_t1)
  MP_pre_t1 = apply(MP_pre_t1,1,prod)


  ### MP PRE T

  for (j in 1:ncol(pre))
  {
    for (i in 1:nrow (pre))
    {
      MP_pre_t[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pre[j],pre[,j][-i],max.pre[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],)}),min.pre[j],max.pre[j]))$value
    }}

  MP_pre_t = 1-(MP_pre_t)
  MP_pre_t = apply(MP_pre_t,1,prod)


  ### Resultado MP PRE
  MP_PRE = (MP_pre_t1)/(MP_pre_t)


  ##############
  ### MP POS ###
  ##############

  ### MP POS T+1

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MP_pos_t1[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pos[j],pos[,j][[i]],max.pos[j],s)}),min.pos[j],max.pos[j]))$value
    }}

  MP_pos_t1 = 1-(MP_pos_t1)
  MP_pos_t1 = apply(MP_pos_t1,1,prod)


  ### MP POS T

  for (j in 1:ncol(pos))
  {
    for (i in 1:nrow (pos))
    {
      MP_pos_t[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min.pos[j],pos[,j][-i],max.pos[j],s))*dpert(x,min.pre[j],pre[,j][[i]],max.pre[j],s)}),min.pre[j],max.pre[j]))$value
    }}

  MP_pos_t = 1-(MP_pos_t)
  MP_pos_t = apply(MP_pos_t,1,prod)


  ### Resultado MC POS
  MP_POS = (MP_pos_t1)/(MP_pos_t)


  #######################
  ###### MP Index #######
  MP = sqrt((MP_PRE)*(MP_POS))


  ###############################
  # PROB. MALMQUIST INDEX (IMP) #
  ###############################

  IMP = sqrt(MC/MP)
  rank = rank(-IMP)
  IMP = cbind(IMP,rank)


  ###############
  ### RESULTS ###
  ###############

  MC = cbind(MC_pre_t,MC_pre_t1,MC_PRE,MC_pos_t,MC_pos_t1,MC_POS,MC)
  MP = cbind(MP_pre_t,MP_pre_t1,MP_PRE,MP_pos_t,MP_pos_t1,MP_POS,MP)

  Result = list(MC=MC, MP=MP, Index=IMP)
  colnames(Result$Index) = c("CPP.Malmquist","Rank")
  Result
}
