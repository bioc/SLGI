\name{gi2005}
\alias{gi2005}
\alias{gi2005.metadata}
\docType{data}
\title{Genetic Interaction Data (EMAP) from the yeast early secretory pathway}
\description{
 The data are in the form of a 424 by 424 array which contains the scores 
 from using the EMAP procedure on yeast strains which are ideally double
 mutants, each strain with a different pair of genes knocked out. For each
 row, the gene named in the row label is knocked out in all pairs, and
 the same holds true for each column.
}
\usage{
 data(gi2005)
 data(gi2005.metadata)
}
\format{
 \code{gi2005} is a 424 by 424 array of real values. 
 \code{gi2005.metadata} is a vector of length 424 which contains
  the common names for the genes that were knocked out. The row and 
 column names of \code{gi2005} are standard names.
}
\details{
 NA values in \code{gi2005} are interactions that were not scored. 
%%FIXME: we might need to find out why, as there is an indication
%%in a related paper that genes that are close in genetic distance are
%%likely to have problems and will always get low scores.
}
\source{
Data were obtained as supplementary material from the publication listed
below.
}
\references{
Schuldiner et al, Exploration of the function and organization of the yeast
early secretory pathway through an epistatic miniarray profile. Cell, 2005,
123:507-519.

Collins et al, A strategy for extracting and analyzing large-scale quantitative
epistatic interaction data. Genome Biology, 2006, 7:R63.
}
\examples{
data(gi2005)
data(gi2005.metadata)
}
\keyword{datasets}
