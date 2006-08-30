\name{predict.pls}
\alias{predict.pls}
\alias{predict.svd}
\title{Classify Observations using Penalized Discriminant Methods }
\description{
  These are functions that can be used to classify new samples (a test
  set) based on an existing classifier created using a training set.
}
\usage{

\method{predict}{pls}(object, x, ...)
\method{predict}{svd}(object, x, ...)
}
\arguments{
  \item{object}{ An object created by a call to \code{pdmClass}. }
  \item{x}{ A matrix of new observations in which rows are samples and
    columns are genes. If not supplied, prediction will be performed on
    the original training set.}
  \item{...}{Other variables passed to predict.}
}

\value{
  A vector of predicted class assignments.
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE}
\author{ Debashis Ghosh}
\examples{
library(fibroEset)
data(fibroEset)
y <- as.numeric(pData(fibroEset)[,2])
x <- t(exprs(fibroEset))
genes <- featureNames(fibroEset)
tmp <- pdmClass(y ~ x)
predict(tmp)
}
\keyword{ models}
\keyword{ robust }
\keyword{ classif}