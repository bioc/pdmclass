\name{pdmGenes}
\alias{pdmGenes}
\title{ A Function to output the Top Ranked Genes from a Penalized
  Discriminant Classifier}
\description{
  After fitting a classifier, it is often desirable to output the most
  "interesting" genes for further validation. This function will output
  the top 'n' genes that discriminate between each class, along with an
  estimate of the stability of the observed rankings (see details for
  more information).
}
\usage{
pdmGenes(formula = formula(data), method = c("pls", "pcr", "ridge"),
data = sys.frame(sys.parent()), weights, theta, dimension = J - 1,
eps = .Machine$double.eps, genelist = NULL, list.length = NULL, B = 100, ...)
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
  \item{dimension}{ The dimension of the solution, no greater than J -
    1, where J is the number of classes. Defaults to J - 1. }
  \item{eps}{ A threshold for excluding small discriminant
    variables. Defaults to \code{.Machine$double.eps}.}
  \item{genelist}{ A vector of gene names, one per gene. }
  \item{list.length}{ The number of 'top' genes to output. }
  \item{B}{ The number of bootstrap samples to use for estimating
    stability. Defaults to 100. More than this may take an inordinate
    amount of time.}
  \item{\dots}{ Additional parameters to pass to \code{method}. }
}
\details{
  The formula interface is identical to all other formula calls in R,
  namely Y ~ X, where Y is a numeric vector of class assignments and X
  is a matrix or data.frame containing the gene expression values. Note
  that unlike most microarray analyses, in this instance the columns of
  X are genes and rows are samples, so most calls will require something
  similar to Y ~ t(X).

  The dimension of the solution is typically J - 1, where J is the
  number of classes. The model fit uses \code{contr.treatment}
  contrasts, which means that all of the coefficients in the model
  are comparing the given class to a baseline class. Therefore, the
  genes listed are those that discriminate between a given class and the
  baseline. For instance, if there are three classes (characterized by a
  numeric vector of 1s, 2s, and 3s), then there will be two sets of 'top
  genes'. The first set will be those genes that discriminate between
  class 2 and class 1, whereas the second set will be the genes that
  discriminate between class 3 and class 1. The 'Y' vector will
  therefore need to be constructed to give the comparisons of interest.
}
\value{
  A list containing a \code{data.frame} for each comparison. The first
  column of each \code{data.frame} contains the gene names, and the
  second column contains the frequency that the gene was observed in the
  bootstrapped samples. 
}
\references{ http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE}
\author{ James W. MacDonald and Debashis Ghosh. Partial least squares
  and principal components regression based on code written by
  Mike Denham and contributed to StatLib. Model fit based on code from
  the \code{mda} package written by Trevor Hastie and Robert Tibshirani
  and ported to R by Kurt Hornik, Brian D. Ripley, and Friedrich Leisch.
}

\examples{
library(fibroEset)
data(fibroEset)
y <- as.factor(pData(fibroEset)[,2])
x <- t(exprs(fibroEset))
genes <- geneNames(fibroEset)
pdmGenes(y ~ x, genelist = genes, list.length = 25, B = 10)
}
\keyword{ models}% at least one, from doc/KEYWORDS
\keyword{ robust }
\keyword{ classif}% __ONLY ONE__ keyword per line
