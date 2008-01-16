\name{seqMatcherAlign}
\alias{seqMatcherAlign}
\alias{getAlignStats}
\title{ Functions to do local alignment of two sequences using EMBOSS matcher }
\description{
  \code{seqMatcherAlign} matches two sequences using the EMBOSS matcher
  program. 

  \code{getAlignStats} extract the statistics from the alignment result data.
}
\usage{
    seqMatcherAlign(pairNameV,BankIDV,seqBank)
    getAlignStats(alignRes)
}
\arguments{
  \item{pairNameV}{ a vector of gene pair names}
  \item{BankIDV}{ a vector of the sequence IDs in the sequence Bank.}
  \item{seqBank}{ a database of all the sequences }
  \item{alignRes}{ object returned by \code{seqMatcherAlign} }
}
\details{
\code{seqMatcherAlign} matches the gene pair names with the sequence
  bank IDs and export the two sequences in to two files: seq1.new and
  seq2.new. Then uses system calls to run EMBOSS matcher program to
  align the two sequences.  The result from matcher is store in file
  "out.matcher".  \code{seqMatcherAlign} read in this file and create a
  R object summarize the alignment results.

  \code{getAlignStats} takes the alignment result data and extract the
  statistics of the result in to data.frame.
}
\value{
  \item{names }{contains the names of the gene pair}
  \item{results }{contains the alignment statistics: the aligned total length,
  the number of identical match, the number of similar match, the number of gaps,
  and the alignment score}
  \item{seq}{displays the aligned sequences}
}
\references{ EMBOSS: The European Molecular Biology Open Software Suite (2000)
  Rice,P. Longden,I. and Bleasby,A.
  Trends in Genetics 16, (6) pp276--277 }
\author{ Z. Jiang }
\note{
  \code{pairMatcherAlign} use system calls to run EMBOSS matcher
  program. You must have EMBOSS matcher installed on your computer.
}

\examples{

seq1 <- "RPHEDEKEAIDEAKKMKVPGENEDESKEEEKSQELEEAIDSKEKSTDARDEQGDEGDNEEENNEEDNENENEHTAPPALVMPSPIEMEEQRM"
seq2 <- "QKYLLKKAIRNFSEYPFYAQNKLIHQQATGLILTEEEKSQELEEKIISKIKKEEHLKKINLKHDYFDLQKKYEKECEILTKLSENLRKEK"
seq3 <- "IHQQATGLILTKIISKIKKEEHVPGENEDLKKINLKHDYFDLQKKYEKECEILTKLSENLRKEEIENKRKEHELMEQKRREEGIETEKEKSLRHPSSSSSSRRRN"
seqBank <- list(seq1=list(seq=seq1),seq2=list(seq=seq2),seq3=list(seq=seq3))
bid <- names(seqBank)
pnames <- c("seq1","seq3")
\dontrun{
ar <- seqMatcherAlign(pnames, bid, seqBank)
ar
getAlignStats(ar)
}
}
\keyword{ methods }
