\name{getInteraction}
\alias{getInteraction}
\title{Count genetic interactions within and between cellular
  organizational units}
\description{
 Count the number of genetic interactions within and between the elements of the interactome. 
}
\usage{getInteraction(iMat, universe, interactome)}
\arguments{
  \item{iMat}{Interaction matrix. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{universe}{Character vector of gene names, e.g., array genes used
    in synthetic genetic array experiments (SGA)}
  \item{interactome}{ Adjacency matrix where row are gene names and
    columns are cellular organizational units. }
}
\value{
  The returned value is a list of 2 matrices:
  \item{bwMat}{A interaction matrix that corresponds to the  cellular
    organizational units interaction matrix where row and columns a
    organizational units names and the value inside the
    matrix are the number of genetic interactions they share.}
  \item{CDs}{Subset of the input interactome that shares interactions.}
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

## Get the interaction
prey <- letters[1:20]
myInteraction <- getInteraction(reducediMat, prey, interactome)
}
\keyword{ methods }
