\name{gi2007}
\alias{gi2007}
\alias{gi2007.metadata}
\docType{data}
\title{ Synthetic Genetic Interaction data from Collins et al}
\description{
 The data \code{gi2007} are a 754 by 754 set of genetic interactions
 that were tested pairwise by either deletion or decreased abundance
 messenger RNA perturbation.
}
\usage{
  data(gi2007)
  data(gi2007.metadata)
}

\format{
  The \code{gi2007} data are a 754 by 754 matrix where values indicate
a score for a synthetic genetic interaction. An \code{NA} indicates that
the genetic interaction was not measured.

  \code{gi2007.metadata} is a data.frame of dimensions 754 rows and
 two columns. The columns are the systematic names
  and the mutation (which is typically either \code{DAMP}, \code{DELETION}
 or the name of the alternate allele that was tested.  In 11 cases an
 alternative allele was tested.
}

\references{
  Collins et al. Nature, 2007, Vol 446, p. 806-810. Data are available
as supplementary material.
}
\examples{
data(gi2007)
data(gi2007.metadata)
}
\keyword{datasets}
