\name{sharedInt}
\alias{sharedInt}

\title{List shared genetic interactions between genes}
\description{
  List shared interactions and cellular organizational units names between genes.
}
\usage{sharedInt(pairL, interactome, threshold=0)}
\arguments{
  \item{pairL}{Dataframe with 3 columns. The first columns are the pair
    of genes tested i.e., the query and array genes. The third columns
    in a logical: TRUE when the 2 genes genetically interact and FALSE
    when they do not.(see \code{AtongPair} dataset as example)}
  \item{interactome}{Adjacency matrix where  row are gene names and
    columns are cellular organizational units names. Each entry has
    value 0 or 1, for absence or presence of a gene in the complex.}
  \item{threshold}{Numeric. Indicate the minimum number of interactions that 2
    genes must share}
}
\value{
 The return value is a list. Each element of the list has for name 2
 genes that genetically interact. Each element of the list corresponds to the
 list of cellular organizational units where the interacting genes are
 found (independently or together).
}

\author{N. LeMeur}
\examples{
## Synthetic genetic interactions 
dat <- data.frame("query" = LETTERS[1:5], "array" = LETTERS[2:6], "interact" = as.logical(sample(c(TRUE, FALSE), 5, TRUE)))
## interactome
interA <- matrix(sample(c(0, 1), 30,TRUE), nrow=6, ncol=5,dimnames = list(LETTERS[1:6], letters[1:5]))

sharedInt(dat, interA, threshold=1)
}
\keyword{data}
\keyword{manip}
