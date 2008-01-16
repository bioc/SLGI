\name{getSharedInteraction}
\alias{getSharedInteraction}
\title{Calculate the number of shared synthetic genetic interactions
  between pairs of genes.}
\description{
 The number of common synthetic genetic
  interacting partners between two genes.
}
\usage{getSharedInteraction(iMat, mode="query")}
\arguments{
  \item{iMat}{Adjacency matrix reporting genetic Interactions. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{mode}{Character vector of value "query" or "target"}
  
}

\value{
  A numeric vector of the number of common genetic interactions between a pair
    of query or target genes.
 }

\author{N. LeMeur}

\seealso{\code{congruence}}

\examples{
intM <- matrix(c(0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0),
                nrow=4, ncol=4, 
                dimnames=list(c("p1","p2","p3","p4"),
                  c("p1","p3","p5","p7")))

sharedInt <- getSharedInteraction(intM)

}
\keyword{ methods }
