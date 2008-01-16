\name{domainDist}
\alias{domainDist}
\title{ Finds the number of gene sets for each shared domain }
\description{
  \code{domainDist} takes a list of shared domains,
  and compute for each distinct domain how many gene sets share it.
}
\usage{domainDist(domainL)}
\arguments{
  \item{domainL}{Each element of the list is a vector of functional
    domains.}
}
\details{
  For each domain that appears in the domain list, \code{domainDist}
  counts the number of elements that have this domain.
}
\value{
  Returns a frequency table with descending order.
}

\author{ Z. Jiang }

\seealso{ \code{\link{getSharedDomains}}, \code{\link{sharedBy}} }
\examples{
data(AtongFnDomain)
domainDist(AtongFnDomain$SharedPfam[1:20])
}
\keyword{ methods }
