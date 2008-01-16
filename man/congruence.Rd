\name{congruence}
\alias{congruence}
\title{Calculate congruence score between pairs of of genes sharing
  pattern of synthetic genetic interactions
  (Ye et al. (2005)).}
\description{
  The \code{congruence} score represents the number of common synthetic genetic
  interacting partners between two genes. The higher is the score
  the more overlap there is between the synthetic genetic partners of
  those genes.
}
\usage{congruence(iMat, sharedInt, mode="query", universe, padjust=FALSE)}
\arguments{
  \item{iMat}{Adjacency matrix reporting genetic Interactions. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column, respectively.}
  \item{sharedInt}{numeric vector representing the number of common
    genetic interactions between a pair
    of query or target genes. See \code{getSharedInteraction}
    for more details}
  \item{mode}{character vector of value "query" or "target"}
   \item{universe}{total number of genes tested}
  \item{padjust}{adjust by the number of genes tested that show at least
  one synthetic genetic interaction.}
}

\value{
  A numeric vector of the congruence score values. 
}
\references{Ye P. et al. (2005). Gene function prediction from congruent
synthetic lethal interactions in yeast. Molecular Systems Biology 1:2005.0026.}

\author{N. LeMeur}

\seealso{\code{getSharedInteraction}}

\examples{
intM <- matrix(c(0,1,0,0,1,1,1,0,1,1,1,1,1,0,1,1),
                nrow=4, ncol=4, 
                dimnames=list(c("p1","p2","p3","p4"),
                  c("p1","p3","p5","p7")))
sharedInt <- getSharedInteraction(intM)
score <- congruence(intM, sharedInt, mode="query", universe=15, padjust=FALSE)
}
\keyword{ methods }
