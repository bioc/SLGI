\name{twoWayTable}
\alias{twoWayTable}
\title{ Generate two-way table for genetic interaction data }
\description{
  Generate two-way table from a vector of genetic interaction status and
  a vector of the pairs that share a functional domain.
}
\usage{
twoWayTable(var1, var2idx)
}
\arguments{
  \item{var1}{ Vector of the status of the first property. }
  \item{var2idx}{ Vector of the index in \code{var1} that have the
    second property. }
}
\details{
  Calculates the count numbers from the given vectors. Then put them
  into a matrix format.
}
\value{
  A two-way contingency table of genetic interaction and whether sharing
  a functional domain.}
\author{ Z. Jiang }
\seealso{ \code{\link{sharedBy}}, \code{\link{getUniquePairs}} }
\examples{
var1 <- c(0,1,1,0,0,0,1,0,1,1)
var2idx <- c(3,5,7)
twoWayTable(var1,var2idx)

data("AtongFnDomain")
pf <- Biobase::reverseSplit(AtongFnDomain$SharedPfam)
idx <- which(rownames(AtongFnDomain$pairs) \%in\% pf$PF00478)
twoWayTable(AtongFnDomain$pairs[,"interact"],idx)
}
\keyword{methods}
