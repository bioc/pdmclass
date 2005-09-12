\name{pdmClass.cv}
\alias{pdmClass.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Leave One Out Crossvalidation }
\description{
  This function performs a leave one out crossvalidation to estimate the
  accuracy of a classifier built using \code{pdmClass}.
}
\usage{
pdmClass.cv(Y, X, method = c("pls", "pcr", "ridge"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ A vector of factors giving the class assignments for the
    samples to be used in the crossvalidation.}
  \item{X}{ A matrix with samples in rows and observations in
    columns. Note that this is different than the usual paradigm for
    microarray data.}
  \item{method}{One of "pls", "pcr", "ridge", corresponding to partial
    least squares, principal components regression and ridge regression.}
}
\details{
  This function performs a leave one out crossvalidation, which can be
  used to estimate the accuracy of a classifier. Each sample is removed
  in turn and a classifier is built using the remaining samples. The
  class of the removed sample is then predicted using the
  classifier. This is repeated for each sample, resulting in a vector of
  predicted class assignments for each sample in the original training
  set.

  Although far from perfect, this method can be used to estimate the
  accuracy of a given classifier without splitting data into a training
  and testing set.
}
\value{
  A vector of factors giving the predicted class assignments for each of
  the samples in the training set. A confusion matrix can be constructed
  using \code{confusion}.
}
\references{http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE

 "Flexible Disriminant Analysis by Optimal Scoring"  by Hastie,
  Tibshirani and Buja, 1994, JASA, 1255-1270.

  "Penalized Discriminant Analysis" by Hastie, Buja and Tibshirani,
  Annals of Statistics, 1995 (in press).}
\author{James W. MacDonald }

\examples{
library(fibroEset)
data(fibroEset)
y <- as.factor(pData(fibroEset)[,2])
x <- t(exprs(fibroEset))
tmp <- pdmClass.cv(y, x)
confusion(tmp, y)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ robust }% __ONLY ONE__ keyword per line
\keyword{classif}