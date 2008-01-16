\name{test2Interact}
\alias{test2Interact}
\title{Summarize genetic interactions within or between cellular
organizational units}
\description{Summarize the genetic interactions within one cellular
organizational unit or between 2 cellular
organizational units. 
}
\usage{test2Interact(iMat, tMat, interactome)}
\arguments{
  \item{iMat}{Genetic interaction matrix. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{tMat}{Adjacency matrix of tested object. Each entry has value 0 or 1,
    representing the fact that the corresponding pairs
    of row and column have been tested for interaction or not.}
  \item{interactome}{Adjacency matrix where  row are gene names and
    columns are cellular organizational units names. Each entry has
    value 0 or 1, for absence or presence of a gene in the complex.}
}
\value{
  the return value is a data.frame with 6 columns.

  \item{unit1, unit2}{cellular organizational units tested and
    interacting}
  \item{tested}{Number of interactions tested between unit1 and unit2}
  \item{interact}{Number of interactions found between unit1 and unit2}
  \item{sizeC1, sizeC2}{Number of genes in  unit1 and unit2}
}
\author{N. LeMeur}

\examples{
set.seed(123)
##Create the interactome
cInt <- sample(c(0,1),30, TRUE)
interactome  <- matrix(cInt, nrow=6, ncol=5,dimnames=list(letters[2:7],LETTERS[1:5]))

## Create cellular organizational units interaction matrix
gInt <- sample(c(1:8), 25, TRUE)
gInt  <- matrix(gInt, nrow=5, ncol=5, dimnames=list(LETTERS[1:5],LETTERS[1:5]))

## All interactome tested
gTest <- matrix(sample(c(0:3), 25, TRUE), nrow=5, ncol=5)
gTested <- gInt+gTest

val <- test2Interact(iMat=gInt, tMat=gTested, interactome=interactome)
}
\keyword{ methods }
