\name{createSquareMatrix}
\alias{createSquareMatrix}


\title{Create a square matrix}
\description{
 Create a square matrix based on row and column names. The new matrix is
 created so that the row and column names are a perfect match and the
 added values are zero.

 In the case of genetic interactions, for example it could be useful
 that the matrix of all the interactions tested and not tested.
}
\usage{createSquareMatrix(data)}

\arguments{
  \item{data}{Matrix}
  }
\value{
  matrix.
}

\author{N. LeMeur}

\examples{
  data(Atong)
  dim(Atong)

  Tong<- createSquareMatrix(Atong)
  dim(Tong)
 }

\keyword{data}
\keyword{manip}