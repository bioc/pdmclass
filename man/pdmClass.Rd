\name{pdmClass}
\alias{pdmClass}
\title{ Function to Classify Microarray Data using Penalized
  Discriminant Methods }
\description{
  This function is used to classify microarray data. Since the
  underlying model fit is based on penalized discriminant methods, there
  is no need for a pre-filtering step to reduce the number of genes.
}
\usage{
pdmClass(formula = formula(data), method = c("pls", "pcr", "ridge"),
data = sys.frame(sys.parent()), weights, theta, dimension = J - 1,
eps = .Machine$double.eps, ...)
}
\arguments{
   \item{formula}{A symbolic description of the model to be fit. Details
    given below. }
  \item{method}{ One of "pls", "pcr", "ridge", corresponding to partial
    least squares, principal components regression and ridge regression.}
  \item{data}{ An optional data.frame that contains the variables in the
    model. If not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{pdmClass} is called. Note that unlike most microarray
    analyses, in this case rows are samples and columns are genes.}
  \item{weights}{An optional vector of sample weights. Defaults to 1. }
  \item{theta}{ An optional matrix of class scores, typically with less
    than J - 1 columns.}
  \item{dimension}{ The dimension of the solution. This will be no
    greater than J - 1 for partial least squares and ridge regression,
    and no greater than J for principal components regression. Defaults
    to J - 1 and J, respectively. }
  \item{eps}{ A threshold for excluding small discriminant
    variables. Defaults to \code{.Machine$double.eps}.}
  \item{\dots}{ Additional parameters to pass to \code{method}. }
}
\details{
  The formula interface is identical to all other formula calls in R,
  namely Y ~ X, where Y is a numeric vector of class assignments and X
  is a matrix or data.frame containing the gene expression values. Note
  that unlike most microarray analyses, in this instance the columns of
  X are genes and rows are samples, so most calls will require something
  similar to Y ~ t(X).
  
\value{
  an object of class \code{"fda"}.  Use \code{predict} to extract
  discriminant variables, posterior probabilities or predicted class
  memberships.  Other extractor functions are \code{coef},
  and \code{plot}.

  The object has the following components:
  \item{percent.explained}{the percent between-group variance explained
    by each dimension (relative to the total explained.)}
  \item{values}{optimal scaling regresssion sum-of-squares for each
    dimension (see reference).  The usual discriminant analysis
    eigenvalues are given by \code{values / (1-values)}, which are used
    to define \code{percent.explained}.}
  \item{means}{class means in the discriminant space.  These are also
    scaled versions of the final theta's or class scores, and can be
    used in a subsequent call to \code{fda} (this only makes sense if
    some columns of theta are omitted---see the references).} 
  \item{theta.mod}{(internal) a class scoring matrix which allows
    \code{predict} to work properly.}
  \item{dimension}{dimension of discriminant space.}
  \item{prior}{class proportions for the training data.}
  \item{fit}{fit object returned by \code{method}.}
  \item{call}{the call that created this object (allowing it to be
    \code{update}-able)}
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE

 "Flexible Disriminant Analysis by Optimal Scoring"  by Hastie,
  Tibshirani and Buja, 1994, JASA, 1255-1270.

  "Penalized Discriminant Analysis" by Hastie, Buja and Tibshirani,
  Annals of Statistics, 1995 (in press).
}
\author{James W. MacDonald and Debashis Ghosh, based on \code{fda} in
  the \code{mda} package of Trevor Hastie and Robert Tibshirani, which
  was ported to R by Kurt Hornik, Brian D. Ripley, and Friedrich Leisch.
}
\examples{
library(fibroEset)
data(fibroEset)
y <- as.numeric(pData(fibroEset)[,2])
x <- t(exprs(fibroEset))
genes <- geneNames(fibroEset)
pdmClass(y ~ x)
}
\keyword{ models}
\keyword{ robust }
\keyword{ classif}