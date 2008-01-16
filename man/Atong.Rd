\name{Atong}
\docType{data}
\alias{Atong}
\alias{tong2004raw}


\title{Systematic genetic analysis with ordered arrays of yeast deletion
  Tong et. al. (2004).}
\description{
  Data from Tong et. al. (2004) buffering experiments using ordered
  arrays of yeast deletion design by Tong et. al. (2001).
}
\usage{
    data(Atong)
    data(tong2004raw)
}
\format{
  \code{tong2004Raw} is dataframe extracted from Table S1 of
   Tong et al. (2004) online supporting material.  We added an extra column,
  queryGene.sysName, which is the  systematic names of the query genes.

  \describe{
    \item{queryGene.geneName}{Column indicates the gene used as query in
      the synthetic genetic array screen (SGA).}
    \item{Int.geneName}{Column indicates the gene identified as an
      interactor with a particular query.}
    \item{Int.sysName}{Column indicates the systematic  name of the open
      reading frame (ORF) that corresponds to the interactor gene.}
    \item{Score}{An interaction scored three times in the three runs by
      visual inspection received a scored of 3. An interaction scored twice
      in the three in the runs by  visual inspection received a scored of 2.
      An interaction scored by the computer-based image analysis but not
      visual inspection received a scored of 1. For interactions that scored
      once in the three runs by visual inspection confirmation was attempted
      only for those genes pairs related functions. Such confirmed
      interactions received a score of 0.}
    \item{RSA}{Column identifies an interaction that was confirmed by random
      spore analysis.}
    \item{Tetrad}{Column identifies an interaction confirmed by tetrad analysis.}
    \item{SS}{Refers to synthetic sick interaction.}
    \item{SL}{Refers to synthetic lethal interaction.}
    \item{Functional.Role}{Column indicates the assigned GO functional
      annotation from their defined subset of annotations.All the
      interactions are identified in this study unless otherwise stated.}
    \item{References}{Genetic Interactions that have been previously
      described.}
    \item{queryGene.sysName}{Column indicates the systematic (ORF) name of
      the gene used as query in a SGA screen.}
  }
  
  \code{Atong} is a 132 by 1008 adjacency matrix of the systematic genetic 
  interactions identified between 132 query genes and the deletion gene
  set (Tong et al. 2001; see \code{\link[SLGI]{SGA}} for more details). 
  The row names correspond to the systematic (ORF) names for the 132 query
  genes. The column names correspond to the systematic (ORF) names of the
  1011 reporter genes, which showed a synthetic lethal or synthetic sick 
  interaction with at least one query genes. Values are 0 or 1, with a 1
  indicating the occurrence of the genetic interaction between the gene pairs.
  }
\source{
  Tong et al, Science. Vol.303, 2004.
}
\references{
  Global Mapping of the Yeast Genetic Interaction Network, Tong et al,
  Science Vol.303, 2004.
  
  \url{http://www.sciencemag.org/cgi/data/303/5659/808/DC1/1}
}
\seealso{
  \code{\link[SLGI]{SGA}}
}

\examples{
 data(Atong)
 dim(Atong)
 }
\keyword{datasets}
