\name{withinComplex}
\alias{withinComplex}

\title{Search for protein co-membership within complexes.}
\description{
  Search for protein co-membership within one (or more) complex(es).
}
\usage{withinComplex(data,interactome)}
\arguments{
  \item{data}{Binary matrix of genes(proteins) linked to other genes(protein) by any
    biological experiment}
  \item{interactome}{Binary matrix composed of genes (rows) and biological
    complexes (columns) \code{\link[ScISI]{ScISI}}}
}

\value{
  Matrix of genes(proteins) co-member of one or more biological complexes.
}

\author{N. LeMeur}

\seealso{\code{\link{byComplex}}}

\examples{
  data(Atong)
  data(ScISIC)
  coMember <- withinComplex(Atong, ScISIC)
  table(coMember)
 }
\keyword{data}
\keyword{manip}
