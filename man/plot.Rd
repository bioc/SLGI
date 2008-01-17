\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,siResult-method}
\alias{plot,siResult,missing-method}

\title{Graphical method to represent the result of the modelSLGI .}
\description{
  a plot method for \code{siResult}.
}
\usage{
\S4method{plot}{siResult}(x,...)
}

\arguments{
  \item{x}{the \code{siResult} object to plot.}
  \item{...}{general commands to be sent to plot.}
}

\details{
  The plot generated from a \code{siResult} object is a dotplot with the
  observed and expected data average of interaction represented in 2
  different colors. 
}

\author{N. LeMeur}

\seealso{\emph{ScISI}}
\examples{
data(ScISIC)
data(Atong)
data(SGA)
model <-  modelSLGI(Atong, universe= SGA, interactome=ScISIC,
type="intM", perm=2)
plot(model)
}
\keyword{data}
\keyword{manip}