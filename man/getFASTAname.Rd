\name{getFASTAname}
\alias{getFASTAname}
\title{ Obtain sequence name from FASTA object }
\description{
  Extract the name of a sequence from a FASTA object that created by
  \code{readFASTA} function from \code{Biostrings} package.
}
\usage{
getFASTAname(Fobj)
}
\arguments{
  \item{Fobj}{ is a FASTA object created by
    \code{\link[Biostrings]{readFASTA}} function from \pkg{Biostrings}}
  package.
} 

\details{
  The function gets the first string between ">" and space the "desc"
  element of the \code{Fobj}, which is the names of the sequence.
}
\value{
  A character string.
}
\author{ Z. Jiang }
\seealso{ \code{\link[Biostrings]{readFASTA}} }
\examples{
f <- gzfile(file.path(.path.package("SLGI"),
            "extdata/orf_trans.fasta.gz"), open = "rt") 
library(Biostrings)
yeastF <- readFASTA(f)
sapply(yeastF[1:5], getFASTAname)
}
\keyword{methods}
