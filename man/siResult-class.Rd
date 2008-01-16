\name{siResult-class}
\docType{class}
\alias{siResult-class}
\alias{show,siResult-method}

\title{A class for representing the result of the SLGI graph permutation
  model.}
\description{A class for representing the result of the \code{modelSLGI} function.}

\section{Slots}{
	 \describe{
	   \item{\code{Observed}:}{Return a \code{"numeric"}
	     vector: the observed  number of synthetic genetic
	 interactions between components of one or two cellular organizational units}
	   \item{\code{Expected}:}{Return a matrix: the expected
	     number of synthetic genetic interactions between
	 components of one or two cellular organizational units}
	 }
}
\section{Methods}{
  \describe{
    \item{\code{plot}}{Graphical representation of the permutation model
      result}
     \item{\code{compare}}{Summarizes the result of the \code{modelSLGI} function}
  }
}
\author{N. LeMeur}

\seealso{
  \code{\link[SLGI]{modelSLGI}},\code{\link[SLGI]{plot}}
}
\examples{
## apply a permutation model 
 data(ScISIC)
 data(Atong)
 data(SGA)
 model <-  modelSLGI(Atong, universe= SGA, interactome=ScISIC,
  type="intM", perm=2)

model

## plot result
plot(model)


}
\keyword{classes}
