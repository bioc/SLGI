\name{normInteraction}
\alias{normInteraction}

\title{Normalize a matrix of biological interactions}
\description{
  Normalize a square matrix of biological interactions according to the number
  of possible interactions between each biological complex.
}
\usage{normInteraction(data, genename, interactome)}
\arguments{
  \item{data}{Square Matrix of biological complexes that shares one or
    more genes(proteins)}
   \item{genename} {Character vector of the gene names that possibly
     create interactions between complexes}
   \item{interactome}{Adjacency matrix where row are genes and
    columns are cellular organizational units. Each entry has
    value 0 or 1, for absence or presence of a gene in a complex, e.g., \code{\link[ScISI]{ScISI}}}
 
}

\value{
  Square matrix of biological complexes linked by one or more
  interacting proteins and normalized by the possible number of
  interactions between each complex.
}

\author{N. LeMeur}

\seealso{\code{\link{getInteraction}}}

\examples{
data(Atong)
data(ScISIC)
data(SGA)
SLa2 <- gi2Interactome(Atong, ScISIC)
## Search for synthetic lethal interaction
compM <- getInteraction(SLa2, SGA, ScISIC)
## Normalize
normIntComplex<- normInteraction(compM$bwMat, SGA, ScISIC)
}
\keyword{data}
\keyword{manip}
