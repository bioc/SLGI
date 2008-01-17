\name{compare}
\docType{methods}
\alias{compare}
\alias{compare,siResult-method}

\title{Compare observed data to expected in permutation models}
\description{
  This method summarizes the result of the \code{modelSLGI} function.
}
\usage{
\S4method{compare}{siResult}(x)
}
\arguments{
  \item{x}{a \code{siResult} object to summarize}
}
\details{
  This compares the number of observed interactions to the number of
expected interactions in each permutation model. It counts how many
times the number of observed interactions is greater than the number of
expected interactions (from the permutations) and divides by the number
of permutations applied.}
\value{
 Numerical vector 
}
\author{N. LeMeur}
\seealso{\emph{modelSLGI}}
\examples{
data(ScISIC)
data(Atong)
data(SGA)
model <-  modelSLGI(Atong, universe= SGA, interactome=ScISIC, type="intM", perm=2)
ans <- compare(model)
}
\keyword{data}
\keyword{manip}
