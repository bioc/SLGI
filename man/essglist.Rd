\name{essglist}
\alias{essglist}
\docType{data}
\title{ The list of yeast essential genes }
\description{
  List of systematic names and common names of the yeast essential genes.
}
\usage{data(essglist)}
\format{
  \code{essglist} is a list with 1103 elements (last download 03/17/2006). The name of each element
  is the systematic gene name. The value of each element is its
  corresponding common (standard) name.  
}
\details{
  The aliases of the yeast gene names can be retrieved with the
  \code{org.Sc.sgdALIAS} environment of the \code{org.Sc.sgd.db} package.
}
\source{
  Saccharomyces Genome Database
  \url{ftp://genome-ftp.stanford.edu/pub/yeast/data_download/literature_curation/phenotypes.tab}
  (last download 03/17/2006)}
\references{
  Saccharomyces Genome Database
  \url{http://www.yeastgenome.org/}
}
\examples{
data(essglist)
essglist[[1]]
names(essglist)
}
\keyword{datasets}
