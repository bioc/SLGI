\name{comemberIn}
\alias{comemberIn}

\title{Retrieve the biological complexes.}
\description{
  Retrieve the biological complexes within which two proteins are comembers.
}
\usage{comemberIn(iMat,interactome)}
\arguments{
  \item{iMat}{Comembership matrix of genes(proteins) that linked to other genes(proteins) by any
    biological experiment}
  \item{interactome}{Adjacency matrix composed of genes (rows) and biological
    complexes (columns) \code{\link[ScISI]{ScISI}}}
}

\value{
  Dataframe of pairs of genes(proteins) and their common biological complexes.
}

\author{N. LeMeur}

\seealso{\code{\link{withinComplex}}}

\examples{
 data(Atong)
 data(ScISI)
 coMember<-withinComplex(Atong,ScISI)
 SLpairWithinComplex <- comemberIn(coMember,ScISI)
 }
\keyword{data}
\keyword{manip}
