\name{AtongPair}
\alias{AtongPair}

\docType{data}
\title{ Data frame the pair of yeast gene tested in Tong et al. 2004. }
\description{
  Data from Tong et. al. buffering experiments (2004) using synthetic
  genetic arrays (SGA)
  (Tong et al. 2001).
}
\usage{data(AtongPair)}
\format{
A data frame with 3 columns and 607881 rows.
}
\details{
  \code{AtongPair} stores the yeast gene names for each tested pairs in Tong
  buffering experiment.  Each row represents one pair.
  \item{query}{Query gene name }
  \item{array}{Array gene name }
  \item{interact}{Logical indicating the synthetic lethal status, if TRUE
    the genetic interaction is lethal.}
}
\source{
  Created from the association matrix reported by Tong et al (2004) and
  the genes from the \code{\link[SLGI]{SGA}} array developed by Tong et
  al. (2001).
}
\references{
 Systematic genetic analysis with ordered arrays of yeast deletion
 mutants, Tong et al, Science. Vol. 294, 2001

  Global Mapping of the Yeast Genetic Interaction Network, Tong et al,
  Science Vol.303, 2004
}

\author{
  N. LeMeur
}

\examples{
data(AtongPair)
dim(AtongPair)
}
\keyword{datasets}
