\name{byComplex}
\alias{byComplex}

\title{Evaluate protein co-membership within cellular organizational units}
\description{
  Count the protein co-members of one (or more) cellular organizational
  units such as complex(es). This co-membership can be characterized by a
  synthetic lethal interaction if bpL is the list of observed synthetic lethal
  interactions or it can be characterized by the number of all the expected interactions within
  that complexes if bpL is all the interactions tested.
}
\usage{byComplex(bpL,interactome)}
\arguments{
  \item{bpL}{List of tested genes (or reported as synthetic lethal) per bait.}
  \item{interactome}{Adjacency matrix where the rows are the genes and the
    columns represent the cellular organizational
  units, e.g., \code{\link[ScISI]{ScISI}}} }

\value{
  Vector of the number of genes(proteins) co-member in one or more
  biological complexes or pathways. 
}

\author{N. LeMeur and R. Gentleman}

\seealso{\code{\link{withinComplex}}}

\examples{
  data(ScISIC)
  data(AtongPair)
  pairSL <-  AtongPair[ AtongPair[,3],]
  SLlist <- split(as.character(pairSL[,2]),as.character(pairSL[,1]))
  ##Number of synthetic lethal pairs within the same complexe 
  bySL <-byComplex(SLlist, ScISIC)
 }
\keyword{data}
\keyword{manip}
