\name{gi2Interactome}
\alias{gi2Interactome}

\title{Reduce genetic interactions matrix}
\description{
  Reduce genetic interactions matrix to the pairs that genetically
  interact and that are present in
  the interactome of interest.
}
\usage{gi2Interactome(iMat, interactome, threshold=0)}
\arguments{
  \item{iMat}{Genetic interaction matrix. Each entry has usually  a
    value of 0 or 1, representing positive or negative interaction of
    corresponding pairs of row and column, respectively.} 
  \item{interactome}{Interactome matrix, e.g. \code{ScISIC}. }
  \item{threshold}{Integer}
}
\value{
  The returned value is the  genetic interaction matrix  reduced to the
  row and column (genes) names that are present in the interactome and where the
  row and column sums are higher than the specified threshold.
}
\author{N. LeMeur}
\examples{
##Create the genetic interaction matrix
gInt <- sample(c(0, 1), 25, TRUE)
iMat  <- matrix(gInt, nrow=5, ncol=5, dimnames=list(letters[1:5],letters[4:8]))

##Create the interactome
cInt <- sample(c(0,1),30, TRUE)
interactome  <- matrix(cInt, nrow=6, ncol=5,dimnames=list(letters[2:7],LETTERS[1:5]))

## Reduce the genetic interaction matrix to match the gene present in
## the interactome
reducediMat <- gi2Interactome(iMat, interactome)
}
\keyword{data}
\keyword{manip}