\name{SGA}
\alias{SGA}
\alias{SGAraw}

\docType{data}
\title{Systematic genetic analysis with ordered arrays of yeast deletion.}
\description{
 Listed of yeast deletion genes used as array probes in the Systematic Genetic
 Analysis (SGA) of yeast deletion Tong et. al. (2001).
}
\usage{
       data(SGAraw)
       data(SGA)
}
\details{

  \code{SGAraw} is a \code{character} vector of length 4672, corresponding to 
  the original yeast deletion genes set on the array. Note that some of those 
  genes correspond to ORFS that have subsequently been rejected.

  \code{SGA} is a \code{character} vector of length 4655, corresponding to 
  the updated yeast deletion genes set on the array. The gene names have
  been updated from common gene name or alias to systemantic names (last
  update Feb. 2006). 

}
\source{
  Table S1 from Tong et al . (2001) online supporting material.
  \url{http://www.sciencemag.org/cgi/content/full/294/5550/2364/DC1}
}
\references{
  Systematic genetic analysis with ordered arrays of yeast deletion
 mutants, Tong et al, Science. Vol. 294, 2001
}
\examples{
 data(SGAraw)
 length(SGAraw)
}
\keyword{datasets}
