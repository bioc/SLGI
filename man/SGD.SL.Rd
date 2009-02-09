\name{SGD.SL}
\alias{SGD.SL}
\alias{SGD.SynRescue}
\alias{SGD.SynGrowthDefect}

\docType{data}
\title{ Interaction data from the Saccharomyces Genome Database}
\description{ The Saccharomyces Genome Database (SGD)
provides, for download a table listing all known interactions in
yeast. This table was downloaded on Jan 25, 2007 and three subsets were
extracted. The synthetic lethal interactions, \code{SGD.SL}, the
synthetic grow defect interactions, \code{SGD.SynGrowthDefect} and the
synthetic rescue interactions, \code{SGD.SynRescue}. No other processing
has been done.  } 

\usage{
  data(SGD.SL)
  data(SGD.SynRescue)
  data(SGD.SynGrowthDefect)
}
\format{
 Each data set is a data frame with  the following 7 variables.
  \describe{
    \item{\code{V1}}{Factor, indicating the type of data.}
    \item{\code{V2}}{Factor describing the interaction,
     in particular naming bait and prey and interactors.}
    \item{\code{V3}}{Factor indicating whether the cells were
     viable.}
    \item{\code{V4}}{Factor which is always NA for these data.}
    \item{\code{V5}}{Factor naming the reference for the interaction.}
    \item{\code{V6}}{Factor with levels indicating the PubMed ID for
     the publication in V5.}
    \item{\code{V7}}{Factor with level \code{BioGRID}, probably
     indicating the source.}
  }
}
\details{
SGD says this about the file:
\preformatted{
Contains interaction data.  Tab-separated columns are:
1) interaction_type (mandatory)
2) genes involved and their mutation type, in the format: ORF
(mutation_type, action), with multiples separated by a | 
3) phenotype (optional, multiples separated by |)
4) description (optional)
5) citation (multiples separated by |)
6) PubMed ID (optional, multiples separated by |)

This file is updated weekly.
}
}
\author{Z. Jiang }
\source{
The file can be downloaded from,
\url{ftp://genome-ftp.stanford.edu/pub/yeast/literature_curation}.
}
\examples{
data(SGD.SL)
}
\keyword{datasets}
