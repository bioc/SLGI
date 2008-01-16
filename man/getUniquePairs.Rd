\name{getUniquePairs}
\alias{getUniquePairs}
\title{ Find unique pairs from an genetic interaction matrix.  }
\description{
  \code{getUniquePairs} can find all the unique pairs from an
  interaction matrix and supplementary array genes, or finds only the
  unique pairs that shows positive interaction.
}
\usage{getUniquePairs(iMat, respV = character(0), only = FALSE)}
\arguments{
  \item{iMat}{Adjacency matrix reporting genetic Interactions. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{respV}{Character vector of all gene names that were tested
    (found to interact or not)}
  \item{only}{ has default value \code{FALSE}, if \code{TRUE}, then only
    reports the positively interacted pairs. }
}
\value{
  A data.frame with two or three columns.  The first two columns are the
  query gene name and the array gene name, respectively.  If \code{only}
  is \code{TRUE}, the third column shows the interaction status. 
}
\author{ Z. Jiang }
\seealso{ \code{\link{getSharedDomains}} }
\examples{
intM <- c(0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0)
dim(intM) <- c(4,4)
dimnames(intM) <- list(c("p1","p2","p3","p4"),c("p1","p3","p5","p7") )
respV <- c("p6","p8")
intM
getUniquePairs(intM,respV,only=FALSE)
getUniquePairs(intM,respV,only=TRUE)
getUniquePairs(intM,only=FALSE)
getUniquePairs(intM,only=TRUE)
}
\keyword{ methods }
