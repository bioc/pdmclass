\name{svdr}
\alias{svdr}
\title{ Classify Microarray Data by Principal Components Regression }
\description{
  This function is used by \code{pdmClass} and \code{pdmGenes} to
  classify microarray data by principal components regression. It is an
  internal function and not intended to be called by the end user.
}
\usage{
svdr(X, y, dimension)
}
\arguments{
  \item{X}{ A numeric matrix assumed to be centred so columns sum to zero.}
  \item{y}{ A numeric vector assumed to sum to zero.}
  \item{dimension}{ The number of PLS factors in the model, which must be
  less than or equal to the rank of \code{X}. }
}
\value{
  An object of class \code{svd}, containing the following items:
  \item{fitted.values }{The fitted values.}
  \item{coefficients }{The model coefficients.}
  \item{dimension}{The number of PLS factors in the model.}
  \item{xmeans}{The sample means.}
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE }
\author{ Debashis Ghosh}
\keyword{internal}

