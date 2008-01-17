\name{modelSLGI}
\alias{modelSLGI}

\title{Permutation model for assessing synthetic genetic interactions in
cellular organizational units.}
\description{
  Permutation model for assessing synthetic genetic interactions within
  and between cellular organizational units such as multi-protein complexes.
}
\usage{
modelSLGI(iMat, universe, interactome,  type="intM", perm=50)
}
\arguments{
  \item{iMat}{Adjacency matrix reporting genetic interactions. Each entry has value 0 or 1,
    representing positive or negative interaction of corresponding pairs
    of row and column.}
  \item{universe}{character vector of the names of the tested
    genes, e.g., names of the genes on the synthetic genetic array
    (SGA) used by Tong et al.}
  \item{interactome}{Adjacency matrix where row are genes and
    columns are cellular organizational units. Each entry has
    value 0 or 1, for absence or presence of a gene in a complex.}
  \item{type}{Character vector of value "intM" (Default) or
    "interactome" to either perform the test based on to the genetic interaction matrix or the
    interactome, respectively.}
  \item{perm}{Number of permutations to apply. Default is 50.}
}

\value{
  Interaction matrix between cellular organizational units.
}

\author{N. LeMeur}

\seealso{\code{\link{getInteraction}}}

\examples{
 data(ScISIC)
 data(Atong)
 data(SGA)
 model <-  modelSLGI(Atong, universe= SGA, interactome=ScISIC,
type="intM", perm=2)
 }
\keyword{data}
\keyword{manip}
