\name{AtongFnDomain}
\alias{AtongFnDomain}
\docType{data}
\title{ The functional domains shared by the tested pairs in Tong et al experiment. }
\description{
  Data developed from Tong et. al. buffering experiments.
}
\usage{data(AtongFnDomain)}
\format{
  A list containing 3 items.
  \item{pairs}{Dataframe of all the gene pairs and their synthetic
    lethality status. }
  \item{SharedPfam}{List of the Pfam domains shared by each pair.
    The order of this list is the same as the order of the pairs.}
  \item{SharedSMART}{List of the SMART domains shared by each pair.
    The order of this list is also the same as the order of the pairs.}
}
\source{
  Created from the association matrix reported by Tong et al. and the
  Pfam (Protein family database \url{http://pfam.janelia.org/}) and SMART database of yeast.
}
\author{Z. Jiang}
\references{
 Global Mapping of the Yeast Genetic Interaction Network, Tong et al,
 Science Vol.303, 2004
 Saccharomyces Genome Database
}
\examples{
data(AtongFnDomain)
names(AtongFnDomain)
}
\keyword{datasets}
