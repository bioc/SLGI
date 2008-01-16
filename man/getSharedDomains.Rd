\name{getSharedDomains}
\alias{getSharedDomains}
\title{ Find domains shared by a given list of gene names. }
\description{
  \code{getSharedDomains} finds domains in the provided environment that
  are shared by a list of genes.
}
\usage{getSharedDomains(geneNameV, env)}
\arguments{
  \item{geneNameV}{ Character vector of gene names.  }
  \item{env}{Name of the domain environment. }
}
\value{
  \code{getSharedDomains} returns a vector of the names of the shared domains.
}
\author{ Z. Jiang }
\seealso{ \code{\link{domainDist}}, \code{\link{sharedBy}} }
\examples{
data(AtongPair)
getSharedDomains(AtongPair["YEL003W-YLR200W",1:2],  YEASTPFAM)
}
\keyword{ methods }
