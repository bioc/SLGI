\name{getTestedPairs}
\alias{getTestedPairs}
\title{ Find interacting and  non-interacting tested pairs from an genetic interaction matrix.}
\description{
  \code{getTestedPairs} find all the pairs from an
  interaction matrix and a list of tested genes.
}
\usage{getTestedPairs(iMat, respV)}
\arguments{
  \item{iMat}{Adjacency matrix reporting genetic Interactions. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{respV}{ Character vector of all gene names that were tested
    (found to interact or not) }
}
\value{
  A data.frame with 4 columns:
  \item{query}{gene names of the query genes}
  \item{array}{gene names og the tested genes (e.g., array genes)}
  \item{interact}{numeric vector of the number of observed interactions
    (0: no interaction; 1: one interaction; 2: two interactions when the
  query genes were also on the array)}
  \item{recip}{logical to indicate whether the reported genes were both
    query and array genes (TRUE: both genes were query and array genes)}. 
}
\author{ N. LeMeur }
\seealso{ \code{\link{getSharedDomains}} \code{\link{getUniquePairs}}}
\examples{
intM <- c(0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0)
dim(intM) <- c(4,4)
dimnames(intM) <- list(c("p1","p2","p3","p4"),c("p1","p3","p5","p7") )
respV <- c("p6","p8")
intM
getTestedPairs(intM,respV)
}
\keyword{ methods }
