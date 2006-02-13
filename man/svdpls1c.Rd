\name{svdpls1c}
\alias{svdpls1c}
\title{Classify Microarray Data by Partial Least Squares }
\description{
  This function is used by \code{pdmClass} and \code{pdmGenes} to
  classify microarray data by partial least squares. It is an
  internal function and not intended to be called by the end user.
}
\usage{
svdpls1c(X, y, ...)
}
\arguments{
  \item{X}{ A numeric matrix assumed to be centred so columns sum to zero.}
  \item{y}{ A numeric vector assumed to sum to zero.}
  \item{\dots}{Further arguments to the function]
  
}
\value{
  An object of class \code{pls}, containing the following items:
  \item{fitted.values }{The fitted values.}
  \item{coefficients }{The model coefficients.}
  \item{dimension}{The number of PLS factors in the model.}
  \item{xmeans}{The sample means.}
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE }
\author{ Mike Denham}
\keyword{internal}