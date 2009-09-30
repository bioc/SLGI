\name{sharedBy}
\alias{sharedBy}
\title{ Find the gene pairs that share a domain. }
\description{
  \code{sharedBy} finds whether the given domain is in each of the
  elements of the domain list.
}
\usage{sharedBy(domainL)}
\arguments{
  \item{domainL}{is a list, each element of the list is a vector of
    domains. }
}
\details{
  \code{sharedBy} first remove all the elements with length 0 or have
  value 'NA'.  Then apply the reverseSplit on the remaining list.
}
\value{
  A list with each element represent a domain, and the values of the
  element are the pairs that share this domain.
}

\author{ Z. Jiang }
\seealso{ \code{\link[Biobase]{reverseSplit}}, \code{\link{domainDist}},
  \code{\link{getSharedDomains}} } 
\examples{
## Load PFAM and SMART domains shared between Tong's Synthetic lethal data
data(AtongFnDomain)
## Find pair that share identical domain
sharedBy(AtongFnDomain$SharedPfam[1:20])
}
\keyword{ methods }
