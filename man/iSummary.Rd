\name{iSummary}
\alias{iSummary}

\title{Summarize cellular organizational units sharing genetic interaction}
\description{
  Summarize the cellular organizational units sharing genetic
  interactions and display their GO annotation if available}
\usage{iSummary(iMat, n=10, reverse=FALSE)}
\arguments{
  \item{iMat}{Comembership matrix of genes(proteins) that linked to
    other genes(proteins) by any biological experiment, e.g., output of
    the \code{getInteraction} function.}
  \item{n}{Numeric threshold indicating the minimum number of genetic
    interactions that a pair of cellular organizational unit must share.}
  \item{reverse}{Logical, by default the function return a list of pair
    of cellular organizational units where the name of each element is
    the number of genetic interactions they share. If reverse is TRUE,
    the output is a vector where the values are the number of
    interactions and the names are the combination of the 2 cellular
    organizational units.}
}
\value{
  The function print the result in the standard output but can also save it
  in variable.
  
  If \code{reverse} is FALSE the output is a list of pairs of cellular
organizational units where the name of each element is the number of
genetic interactions they share.
  
If \code{reverse} is TRUE the output is a vector where the values are
the number of interactions and the names are the combination of the 2
cellular organizational units.  }

\author{N. LeMeur}
\examples{
data(Atong)
data(ScISIC)
data(SGA)
SLa2 <- gi2Interactome(Atong, ScISIC)
## Search for synthetic lethal interaction
compM <- getInteraction(SLa2, SGA, ScISIC)
## Display the tightly interacting pairs
largeInt <- iSummary(compM$bwMat,n=15)
}
\keyword{data}
\keyword{manip}
