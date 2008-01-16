\name{dSLAM.GPL1444}
\alias{dSLAM.GPL1444}
\alias{dSLAM}
\docType{data}
\title{dSLAM platform used for Synthetic Lethal screens in the Boeke Lab}
\description{
 These data are the 21991 probes spotted on the dSLAM array (heterozygote
 diploid-based synthetic lethality analyzed by microarray) used to test
 synthetic lethal interactions by Pan et al (2006).
}
\usage{
data(dSLAM.GPL1444)

data(dSLAM)
}
\format{
  dSLAM.GPL1444 is a data frame with 21991 observations on the following 10 variables.
  \describe{
    \item{\code{ID}}{Serial identifier for probe.}
    \item{\code{ROW}}{Row number in the array as scanned with GenePix scanner.}
    \item{\code{COLUMN}}{Column number in the array as scanned with GenePix scanner.}
    \item{\code{TAGTYPE}}{Code for whether tag is 5' (Up) or 3' (Dn) relative to the open reading frame (ORF).}
    \item{\code{PROBE}}{Code for singleton probes arrayed in ORF order (ArrA, ArrB), five-fold replicate probes arrayed in randomized order (Rpts), systematic mutations arrayed across the center of the array (Muts), negative controls (NegT), or probes peripheral to the array as specified by the manufacturer (Edge)}
    \item{\code{ORF}}{Systematic ORF name (from SGD, Feb 2003) (\url{http://genome-www4.stanford.edu/cgi-bin/SGD/locus.pl?locus=}}
    \item{\code{GENE}}{Standard gene name (SGD) (or ORF if not available)}
    \item{\code{SEQUENCE}}{DNA sequence of probe (includes custom-designed sequences for 193 YA* and YM* ORFs missing DnTags)}
    \item{\code{SGDID}}{Unique ORF identifier from SGD; 'S000000000' denotes missing value}
    \item{\code{SPOT_ID}}{spot identifier; ('YQL' ORFs denote custom-designed sequences; 'NegA', 'NegB', 'PosA', 'PosB' denote proprietary sequences specified by the manufacturer)}
  }

  dSLAM is a character vector of length 5641 that contains the unique and valid systematic ORF names. 
}
\details{
  The dSLAM.GPL1444 were directly obtain from parsing the
  GPL1444\_family.soft.gz available at
\url{http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GPL1444}

  dSLAM is a vector of length 5641, extracted from the dSLAM.GPL1444
  ORF, and that contains the unique and valid systematic ORF names. This
  vector was built in three steps. First the ORFs with SGDID equals to
  S000000000 in the dSLAM.GPL1444 data frame were removed as some
  correspond to custom sequences and other were dubious ORFs that have
  been  deleted from SGD or merged with other ORFs. Secondly, the duplicated
  names were removed.  Then, the systematic ORF names were verified against
  the YEAST data package.
  
}
\source{
The data were extracted from the Gene Expression Omnibus (GEO) website: 
  \url{http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GPL1444}
}
\references{
Pan X, Ye P, Yuan DS, Wang X, Bader JS, Boeke JD. 
A DNA integrity network in the yeast Saccharomyces cerevisiae.
Cell. 2006 Mar 10;124(5):1069-81
}
\seealso{\code{\link[SLGI]{Boeke2006raw}} \code{\link[SLGI]{Boeke2006}}}

\keyword{datasets}
