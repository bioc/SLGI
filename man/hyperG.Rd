\name{hyperG}
\alias{hyperG}

\title{Hypergeometric test}
\description{
  a hypergeometric test for genetic interaction data.
}
\usage{hyperG(data, nbTested, universe)}
\arguments{
  \item{data}{Matrix with 2 columns the first one corresponds to the
    number of interactions per pair of interacting complexes
    and the second one to number of tested interactions. This could be the first
    two columns resulting from a call to the \code{test2Interact} function.}
  \item{nbTested}{Number of interacting pairs}
  \item{universe}{Total Number of tested pairs}
}

\seealso{
  \code{phyper}
}
\author{N. LeMeur}
\examples{
## Create matrix interaction x tested matrix
interact <-  c(1, 3, 2, 2, 6, 5, 2, 4, 1, 3)
tested <-  c(3, 3, 5, 4, 8, 5, 3, 4, 2, 3)
mat <- cbind(interact, tested)

## Perform test
res <- hyperG(mat, 1000, 10000)
summary(res$P)
}
\keyword{data}
\keyword{manip}