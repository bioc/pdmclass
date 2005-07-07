\name{my.svd}
\alias{my.svd}
\title{A Function to Perform Singular Value Decomposition }
\description{
  An alternative to Singular Value Decomposition function svd
  that examines n by p matrix x and if n < p obtains the svd 
  by applying svd to the transpose of x. This is an internal function
  and is not intended to be called by the end user.
}
\usage{
my.svd(x, nu = min(n, p), nv = min(n, p))
}
\arguments{
  \item{x}{ A numeric or complex matrix }
  \item{nu}{ The number of left singular vectors to be computed.}
  \item{nv}{ The number of right singular vectors to be computed. }
}
\details{
 This implementation of SVD uses the LINPACK routines DSVDC for numeric
 matrices and ZSVDC for complex matrices.
}
\value{
  The returned value is a list with components:
  \item{d }{A vector containing the singular values of \code{x}}
  \item{u }{A matrix whose columns contain the left singular vectors of
    \code{x}, present if 'nu > 0'.}
  \item{v}{A matrix whose columns contain the right singular vectors of
    \cose{x}, present if 'nv > 0'.}
}
\references{  http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE}
\author{Mike Denham}

\keyword{ internal}
