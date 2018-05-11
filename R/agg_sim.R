#' Aggregation of expert's estimatives by similarity of values
#' @description This function computes the aggregated value of different expert's estimatives, using Beta PERT distributions to randomize the decision matrix.
#' @param x Decision matrix of expert estimatives (rows) and criteria (columns). Benefit criteria must be positive and cost criteria must be negative.
#' @param min Vector of minimum values in each criterion scale. For common scales to all criteria, the vector must repeat the minimum value as many times as the number of criteria.
#' @param max Vector of maximum values in each criterion scale. For common scales to all criteria, the vector must repeat the maximum value as many times as the number of criteria.
#' @param s Shape of a Beta PERT distribution, as described in the package 'mc2d'. There is no default value, however the higher the shape the higher the kurtosis of the random variable.
#' @param w Weights describing the expert experience in the subject matter.
#' @param b Beta describes the balance between the expert weights and their opinions. Beta varies in the interval [0,1]. The higher the index, the higher the importance of weights.
#' @return SM are the Similarity Matrices per criterion. CDC describes the Consensus Coefficient matrix. Agg.value gives the aggregated value of expert opinions per criterion.
#' @examples
#' ## Expert's estimatives on four criteria
#' Exp.1 = c(4,7,6,8)
#' Exp.2 = c(4,3,6,5)
#' Exp.3 = c(3,8,2,9)
#' Exp.4 = c(6,8,9,7)
#' Exp.5 = c(5,9,2,4)
#' Exp.6 = c(7,6,5,5)
#' x = rbind(Exp.1,Exp.2,Exp.3,Exp.4,Exp.5) # Decision matrix
#' min = c(0,0,0,0) # Minimum scale values.
#' max = c(10,10,10,10) # Maximum scale values.
#' s = 4 # Shape
#' w = c(0.4,0.3,0.2,0.06,0.04) # Expert relevance.
#' b = 0.4
#' Agg.Sim(x,min,max,s,w,b)
#' @importFrom mc2d dpert ppert
#' @importFrom stats integrate
#' @export
Agg.Sim = function(x,min,max,s,w,b){
PMax = x
A = x
B = x

for (j in 1:ncol(x))
{
  for (i in 1:nrow(x))
  {
    PMax[i,j] = (integrate(Vectorize(function(x) {prod(ppert(x,min[j],A[,j][-i],max[j],s))*dpert(x,min[j],B[,j][[i]],max[j],s)}),min[j],max[j])) $value
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
    PMin[i,j] = (integrate(Vectorize(function(x) {prod(1-ppert(x,min[j],A[,j][-i],max[j],s))*dpert(x,min[j],B[,j][[i]],max[j],s)}),min[j],max[j])) $value
  }}
PMin = PMin[,]

# Number of Experts for looping
NrEsp = nrow(x)
A = x
B = x

# Cross matrices PMax e Pmin

PMax = vector("list", ncol(x))
PMin = vector("list", ncol(x))

### P.MAX e P.min

for (k in 1:(ncol(x)))
{
  PMax[[k]] = matrix(0,nrow=NrEsp,ncol=NrEsp)
  PMin[[k]] = matrix(0,nrow=NrEsp,ncol=NrEsp)
  for (j in 1:(NrEsp))
  {
    for (i in 1:(NrEsp))
    {
      PMax[[k]][i,j] = (integrate(Vectorize(function(x) {(ppert(x,min[k],A[i,k],max[k],s))*dpert(x,min[k],B[j,k],max[k],s)}),min[k],max[k]))$value
      PMin[[k]][i,j] = (integrate(Vectorize(function(x) {(1-ppert(x,min[k],A[i,k],max[k],s))*dpert(x,min[k],B[j,k],max[k],s)}),min[k],max[k]))$value
    }}}

# Similarity matrices per criterion

S = vector("list", ncol(x))
for (k in 1:(ncol(x)))
{
  S[[k]] = PMax[[k]]/PMin[[k]]
  S[[k]] = ifelse(S[[k]]>1,1/S[[k]],S[[k]])
  S[[k]] = ifelse(S[[k]]>0.999,1,S[[k]])
  rownames(S[[k]]) = paste0("Exp",1:NrEsp)
  colnames(S[[k]]) = paste0("Exp",1:NrEsp)
}
names(S) = paste0("Crit",1:ncol(x))

### Aggregation of estimatives

AE = vector("list", ncol(x))
for (k in 1:(ncol(x)))
{
  AE[[k]] = apply(S[[k]],1,sum)
  AE[[k]] = (AE[[k]]-1)/(NrEsp-1)
}
names(AE) = paste0("Crit",1:ncol(x))

### Relative Aggregation (RAD) and Coeficient of Consensus (CDC)

Soma = vector("list", ncol(x))
RAD = vector("list", ncol(x))
CDC = vector("list", ncol(x))
for (k in 1:(ncol(x)))
{
  Soma[[k]] = sum(AE[[k]])
  RAD[[k]] = AE[[k]]/Soma[[k]]
  CDC[[k]] = b*w+(1-b)*RAD[[k]]
}
names(CDC) = paste0("Crit",1:ncol(x))
CDC = t(matrix(unlist(CDC),ncol(x),NrEsp))
rownames(CDC) = paste0("Exp",1:NrEsp)
colnames(CDC) = paste0("Crit",1:ncol(x))

### Final aggregation (New alternative parameters)

A.mode = vector("list", ncol(x))
for (k in 1:(ncol(x)))
{
  A.mode[[k]] = sum(CDC[[k]]*A[,k])
}
names(A.mode) = paste0("Crit",1:ncol(x))
A.mode = unlist(A.mode)

Result = list(SM = S, CDC = CDC, Agg.value = A.mode)
Result
}
