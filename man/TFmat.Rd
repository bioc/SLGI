\name{TFmat}
\alias{TFmat}
\docType{data}
\title{ Transcription Factor Binding Affinities}
\description{
The data are from Lee et al, the rows of the matrix represent genes in
S. cerevisiae, the columns known transcription factor. The value in each
entry represents the p-value, as reported by Lee et al, for the
transcription factor (TF) binding upstream of the gene.
}
\usage{data(TFmat)}
\format{
\code{TFmat} is a matrix, rows represent genes, columns transcription
factors and the elements are p-values representing some notion of the
likelihood that the transcription factor binds up stream of the gene.
}
\author{Z. Jiang}
\source{
  Supplementary material from
  \url{http://web.wi.mit.edu/young/regulator_network/}
}
\references{
  Transcriptional Regulatory Networks in Saccharomyces cerevisiae,
  Lee et al. Science 298:799-804 (2002).
}
\examples{
data(TFmat)
}
\keyword{datasets}
