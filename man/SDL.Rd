\name{SDL}
\alias{SDL}
\alias{SLchr}

\docType{data}
\title{ The Association matrix for the synthetic dosage
lethal screens in Yeast. }
\description{
The data reported in Table 6 of the supplementary data of
Measday et. al.
}
\usage{
       data(SDL)
       data(SLchr)
}
\format{
\code{SDL} is a matrix with 141 rows and 9 columns. The columns represent 
3 genes at each of 3 temperatures (16, 25, 37 Celsius). The gene names and temperatures are 
combined in the column names. The row names are yeast standard names.
The values are \code{NA}, no effect, \code{SDS} for synthetic dosage sick,
\code{SL} for synthetic lethal and \code{SDL} for synthetic dosage lethal.

\code{SLchr} is a matrix with 84 rows and 14 columns. Each column represents
a query strain which was tested against the genome wide set of deletion
strains. The entries can be \code{NA} for no effect, \code{SL} for synthetic
lethal and \code{SS} for synthetic sick.

}

\source{
Supplementary Table 6 of the reference given below.
}
\references{
Systematic yeast synthetic lethal and synthetic dosage lethal screens identify
genes required for chromosome segregation. Measday et al, PNAS, 2005,
13956-13961.
}
\examples{
data(SDL)
table(SDL)
}
\keyword{datasets}
