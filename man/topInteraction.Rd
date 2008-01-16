\name{topInteraction}
\alias{topInteraction}

\title{Extract interacting biological complexes}
\description{
  Extract the top X interacting biological complexes.
}
\usage{topInteraction(data,top=10)}
\arguments{
  \item{data}{Square matrix of biological complexes that shares one or
    more genes(proteins)}
  \item{top}{Interger that represents the percentage of interacting complexe}

}

\value{
  Data frame of biological complexes that interact. The first two
  columns are the cellular organizational units names and the third
  column indicates the number of interactions.
}

\author{N. LeMeur}

\examples{
data(Atong)
data(ScISIC)
data(SGA)
SLa2 <- gi2Interactome(Atong, ScISIC)
## Search for synthetic lethal interaction
compM <- getInteraction(SLa2, SGA, ScISIC)
top10Interaction<- topInteraction(compM$bwMat,top=10)
}
\keyword{data}
\keyword{manip}
