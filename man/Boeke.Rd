\name{Boeke}
\alias{Boeke2006raw}
\alias{Boeke2006}
\docType{data}
\title{Incidence matrix of Synthetic Lethal interaction from the Boeke Lab}
\description{
 Data from Pan et al. experiments on DNA integrity Network in the
 Yeast S. cerevisiae.
}
\usage{
data(Boeke2006raw)
data(Boeke2006)
}
\format{
  \code{Boeke2006raw} is a data frame with 5775 observations on the following 6 variables.
  \describe{
    \item{\code{Query.ORF}}{ORF associated with the query gene.}
    \item{\code{Query.Gene}}{Common name of the query gene.}
    \item{\code{Target.ORF}}{ORF for the array gene.}
    \item{\code{Target.Gene}}{Common name of the array gene.}
    \item{\code{RSA}}{Random spore analysis}
    \item{\code{Tetrad}}{Tetrad dissection}
  }

  \code{Boeke2006} is an incidence matrix is a 74 by 843 adjacency matrix of the systematic genetic 
  interactions identified between 74 query genes and the deletion gene
  set in Pan et al.(2004). The row names correspond to the systematic
  (ORF) names for the 74 query genes. The column names correspond to the
  systematic (ORF) names of the
  843 reporter genes, which showed a synthetic lethal or synthetic sick 
  interaction with at least one query genes. Values are 0 or 1, with a 1
  indicating the occurrence of the genetic interaction between the gene
  pairs.
}
\details{
In Pan et al (2006), the authors provide this note. Note: SL - synthetically lethal; 
SF/SL-very severe synthetic fitness defects; 
SF-obvious but modest synthetic fitness defects; 
SF (slight) - slight synthetic fitness defect. Approximately 10\% of the 
positive interactions presented here were not scored as positive in the 
dSLAM screens.  These were individually tested because we wanted to make 
sure that they were indeed false negatives in the dSLAM screens.  We also 
note that there is a small chance that the interactions scored as positive in 
RSA (random spore analysis) might not reflect direct growth defects of the 
double mutants but rather, the double mutants are defective in expressing the 
MFA1pr-HIS3 reporter. 
}
\source{
The data were extracted from Pan et al (2004) and Table S1 of Pan et al. (2006).
}
\references{
Pan X, Ye P, Yuan DS, Wang X, Bader JS, Boeke JD. 
A DNA integrity network in the yeast Saccharomyces cerevisiae.
Cell. 2006 Mar 10;124(5):1069-81

Pan X, Yuan DS, Xiang D, Wang X, Sookhai-Mahadeo S, Bader JS, Hieter P,
Spencer F, Boeke JD.
A robust toolkit for functional profiling of the yeast genome.
Mol Cell. 2004 Nov 5;16(3):487-96.
}
\seealso{\code{\link[SLGI]{dSLAM.GPL1444}}, and \code{\link{dSLAM}}}


\keyword{datasets}
