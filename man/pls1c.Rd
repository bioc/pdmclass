\name{pls1c}
\alias{pls1c}
\title{ Function to Fit Modified Helland Algorithm}
\description{
  This function fits a partial least squares model based on the modified
  Helland algorithm. This is an internal function and is not intended to
  be called by the end user.
}
\usage{
pls1c(X, y, dimension = min(dx[1] - 1, dx[2]))
}
\arguments{
  \item{X}{ A numeric matrix assumed to have been centred so columns sum
  to zero.}
  \item{y}{ A numeric vector assumed to sum to zero.}
  \item{dimension}{The number of PLS factors in the model, which must be
  less than or equal to the rank of \code{X}.}
}
\value{
  An object of class \code{pls}, containing the following items:
  \item{fitted.values }{The fitted values.}
  \item{coefficients }{The model coefficients.}
  \item{dimension}{The number of PLS factors in the model.}
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE }
\author{ Mike Denham }
\keyword{ internal }