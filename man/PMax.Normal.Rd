% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmax_normal.R
\name{PMax.Normal}
\alias{PMax.Normal}
\title{Probabilities of maximization, by Normal distributions}
\usage{
PMax.Normal(x)
}
\arguments{
\item{x}{Decision matrix of Alternatives (rows) and Criteria (columns). Benefit criteria must be positive and cost criteria must be negative.}
}
\value{
PMax are the joint probabilities of each alternative being higher than the others, per criterion.
}
\description{
This function computes the probabilities of each alternative maximizing the preference per criterion, using Normal distributions to randomize the decision matrix.
}
\examples{
# Decision matrix
Alt.1 = c(2,30,86,-5)
Alt.2 = c(4,26,77,-12)
Alt.3 = c(3,22,93,-4)
Alt.4 = c(6,34,65,-10)
Alt.5 = c(5,31,80,-8)
x = rbind(Alt.1,Alt.2,Alt.3,Alt.4,Alt.5)
PMax.Normal(x)
}
\references{
Sant'Anna, Annibal P. (2015). Probabilistic Composition of Preferences: Theory and Applications, Springer.
}
