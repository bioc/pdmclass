##########################################################
##
##
## pdmclass - A package for classification of microarray data
##            based on penalized discriminant methods.
##
##  Ghosh, D., Penalized discriminant methods for the classification
##  of tumors from gene expression data. Biometrics, 59, 992-1000(2003).
##
##
## History:
## 6-28-05 First pass with functions from http://www.sph.umich.edu/~ghoshd/COMPBIO/POPTSCORE/POSmacros.txt
## 7-6-05  Added functionality to output top gene list
##
##############################################################


"pls1c"<-
  function(X, y, dimension = min(dx[1] - 1, dx[2]))
{
  ## Copyright Mike Denham, October 1994.
  ## Comments and Complaints to: snsdenhm@reading.ac.uk
  ##
  ## Modified Helland Algorithm (Helland 1988 + Denham 1994)
  ##
  ## X: A matrix which is assumed to have been centred so that columns
  ##    sum to zero.
  ##
  ## y: A vector assumed to sum to zero.
  ##
  ## dimension: The number of PLS factors in the model which must be less than or
  ##    equal to the  rank of X.
  ##
  ## Returned Value is the vector of PLS regression coefficients
  ##
  ## tol is set as the tolerance for the QR decomposition in determining
  ## rank deficiency
  ##
  tol <- 1e-10
  X <- as.matrix(X)
  dx <- dim(X)
  W <- matrix(0, dx[2], dimension)
  XW <- matrix(0, dx[1], dimension)
  s <- crossprod(X, y)
  W[, 1] <- s
  XW[, 1] <- X %*% s
  QR <- qr(XW[, 1], tol = tol)
  r <- qr.resid(QR, y)
  if(dimension > 1) {
    for(i in 2:dimension) {
      w <- crossprod(X, r)
      W[, i] <- w
      XW[, i] <- X %*% w
      QR <- qr(XW[, 1:i], tol = tol)
      r <- qr.resid(QR, y)
    }
  }
  coef <- W %*% qr.coef(QR, y)
  fitted <- X %*% coef
  structure(list(fitted.values = fitted, coefficients = coef, dimension = dimension),
            class="pls") 
}

"svdpls1c"<-
  function(X, y, dimension = r)
{
  ## Copyright Mike Denham, October 1994.
  ## Comments and Complaints to: snsdenhm@reading.ac.uk
  ##
  ## Modified Helland Algorithm (Helland 1988 + Denham 1994)
  ##
  ## X: A matrix which is assumed to have been centred so that columns
  ##    sum to zero.
  ##
  ## y: A vector assumed to sum to zero.
  ##
  ## dimension: The number of PLS factors in the model which must be less than or
  ##    equal to the  rank of X.
  ##
  ## Returned Value is the vector of PLS regression coefficients
  ##
  ## tol is set as the tolerance for the QR decomposition in determining
  ## rank deficiency
  ##

  tX <- as.matrix(X)
  r <- min(dim(X) - c(1, 0))
  X <- svd(X)
  fitted <- NULL
  coef <- NULL
  mm <- colMeans(tX)
  dy <- dim(y)[2]
  for (i in 1:dy) {
    tmp <- pls1c(diag(X$d[1:r]), crossprod(X$u[, 1:r], y[,i]), dimension)
    tcoef <- X$v[, 1:r] %*% tmp$coefficients
    coef <- cbind(coef,tcoef)
    fitted <- cbind(fitted,tX %*% tcoef)
  }
  structure(list(fitted.values=fitted,coefficients=coef,xmeans=mm,dimension=dimension),
            class="pls")
}

"svdr" <- function(X, y, dimension) {

  r  <- dimension
  tX <- as.matrix(X)
  tsvd <- my.svd(X)
  mm <- apply(X,2,mean)
  fitted <- NULL
  coef <- NULL
  dy <- dim(y)[2]
  for (i in 1:dy) {
    if (r == 1) 
      tmp <- lm.fit(as.matrix(tsvd$u[,1:r]) %*% tsvd$d[1],y[,i])
    else
      tmp <- lm.fit(as.matrix(tsvd$u[,1:r]) %*% as.matrix(diag(tsvd$d[1:r])),y[,i])
    tcoef <- as.matrix(tsvd$v[,1:r]) %*% tmp$coefficients
    coef <- cbind(coef, tcoef)
    fitted <- cbind(fitted, tX %*% tcoef)
  }			
  structure(list(fitted.values = fitted, coefficients = coef, dimension = dimension,
                 xmeans = mm), class="svd")

}


"my.svd"<-
  function(x, nu = min(n, p), nv = min(n, p))
{
  ## Alternative to Singular Value Decomposition function svd
  ## Examines matrix n by p matrix x and if n < p obtains the svd 
  ## by applying svd the transpose of x.
  x <- as.matrix(x)
  dmx <- dim(x)
  n <- dmx[1]
  p <- dmx[2]
  transpose.x <- n < p
  if(transpose.x) {
    x <- t(x)
    hold <- nu
    nu <- nv
    nv <- hold
  }
  cmplx <- mode(x) == "complex"
  if(!(is.numeric(x) || cmplx))
    stop("x must be numeric or complex")
  if(!cmplx)
    storage.mode(x) <- "double"
  dmx <- dim(x)
  n <- dmx[1]
  p <- dmx[2]
  mm <- min(n + 1, p)
  mn <- min(dmx)
  job <- (if(nv) 1 else 0) + 10 * (if(nu == 0) 0 else if(nu == mn)
                                   2
  else if(nu == n)
                                   1
  else stop("Invalid value for nu (must be 0, number of rows, or number of cols)"
            ))
  z <- .Fortran(if(!cmplx) "dsvdc" else "zsvdc",
		x,
		as.integer(n),
		as.integer(n),
		as.integer(p),
		d = if(!cmplx) double(mm) else complex(mm),
		if(!cmplx) double(p) else complex(p),
		u = if(!cmplx) if(nu)
                matrix(0, n, nu)
                else 0 else if(nu)
                matrix(as.complex(0), n, nu)
		else as.complex(0),
		as.integer(n),
		v = if(!cmplx) if(nv)
                matrix(0, p, p)
                else 0 else if(nv)
                matrix(as.complex(0), p, p)
		else as.complex(0),
		as.integer(p),
		if(!cmplx) double(n) else complex(n),
		as.integer(job),
		info = integer(1),
                PACKAGE = "base")[c("d", "u", "v", "info")]
  if(z$info)
    stop(paste("Numerical error (code", z$info, ") in algorithm"))
  if(cmplx) {
    if(all(Im(z$d) == 0))
      z$d <- Re(z$d)
    else stop("a singular value has a nonzero imaginary part")
  }
  length(z$d) <- mn
  if(nv && nv < p)
    z$v <- z$v[, seq(nv)]
  if(transpose.x) {
    z <- z[c("d", if(nu) "u" else NULL, if(nv) "v" else NULL)]
    names(z) <- names(z)[c(1, 3, 2)]
    z
  }
  else {
    z[c("d", if(nv) "v" else NULL, if(nu) "u" else NULL)]
  }
}

predict.pls <- function(object, x, ...) {

    if (missing(x)) 
        fitted(object)
    else scale(x, object$xmeans, FALSE) %*% object$coef
}

predict.svd <- function(object, x, ...) {

    if (missing(x)) 
        fitted(object)
    else scale(x, object$xmeans, FALSE) %*% object$coef
}


pdmClass <- function (formula = formula(data), method = c("pls", "pcr", "ridge"),
                      data = sys.frame(sys.parent()), weights, theta,
                      dimension = J - 1, eps = .Machine$double.eps, ...){
  
  this.call <- match.call()
  m <- match.call(expand = FALSE)
  m[[1]] <- as.name("model.frame")
  m <- m[match(names(m), c("", "formula", "data", "weights"), 
               0)]
  m <- eval(m, sys.frame(sys.parent()))
  Terms <- attr(m, "terms")
  g <- model.extract(m, response)
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  dd <- dim(x)
  n <- dd[1]
  weights <- model.extract(m, weights)
  if (!length(weights)) 
    weights <- rep(1, n)
  else if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(g) != n) 
    stop("g should have length nrow(x)")
  fg <- factor(g)
  prior <- table(fg)
  prior <- prior/sum(prior)
  cnames <- levels(fg)
  g <- as.numeric(fg)
  J <- length(cnames)
  iswt <- FALSE
  if (missing(weights)) 
    dp <- table(g)/n
  else {
    weights <- (n * weights)/sum(weights)
    dp <- tapply(weights, g, sum)/n
    iswt <- TRUE
  } 
  if (missing(theta)) 
    theta <- contr.helmert(J)
  theta <- contr.fda(dp, theta)
  Theta <- theta[g, , drop = FALSE]
  method <- match.arg(method)
  fit <- switch(method,
                pls = svdpls1c(scale(x, center = TRUE, scale = FALSE),
                  scale(Theta, center = TRUE, scale = FALSE), dimension),
                pcr = svdr(x, Theta, dimension + 1),
                ridge = gen.ridge(x, Theta, weights, lambda = 1, ...))
  if (iswt) 
    Theta <- Theta * weights
  ssm <- t(Theta) %*% fitted(fit)/n
  ed <- svd(ssm, nu = 0)
  thetan <- ed$v
  lambda <- ed$d
  lambda[lambda > 1 - eps] <- 1 - eps
  discr.eigen <- lambda/(1 - lambda)
  pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
  dimension <- min(dimension, sum(lambda > eps))
  if (dimension == 0) {
    warning("degenerate problem; no discrimination")
    return(structure(list(dimension = 0, fit = fit, call = this.call), 
                     class = "fda"))
  }
  thetan <- thetan[, seq(dimension), drop = FALSE]
  pe <- pe[seq(dimension)]
  alpha <- sqrt(lambda[seq(dimension)])
  sqima <- sqrt(1 - lambda[seq(dimension)])
  vnames <- paste("v", seq(dimension), sep = "")
  means <- scale(theta %*% thetan, FALSE, sqima/alpha)
  dimnames(means) <- list(cnames, vnames)
  names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
  names(pe) <- vnames
  obj <- structure(list(percent.explained = pe, values = lambda, 
                        means = means, theta.mod = thetan, dimension = dimension, 
                        prior = prior, fit = fit, call = this.call, terms = Terms), 
                   class = "fda")
  obj
}



pdmGenes <- function (formula = formula(data), method = c("pls", "pcr", "ridge"),
                      data = sys.frame(sys.parent()), weights, theta,
                      dimension = J - 1, eps = .Machine$double.eps,
                      genelist = NULL, list.length = NULL,
                      B = 100, ...){
  
  this.call <- match.call()
  m <- match.call(expand = FALSE)
  m[[1]] <- as.name("model.frame")
  m <- m[match(names(m), c("", "formula", "data", "weights"), 
               0)]
  m <- eval(m, sys.frame(sys.parent()))
  Terms <- attr(m, "terms")
  g <- model.extract(m, response)
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  dd <- dim(x)
  n <- dd[1]
  weights <- model.extract(m, weights)
  if (!length(weights)) 
    weights <- rep(1, n)
  else if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(g) != n) 
    stop("g should have length nrow(x)")
  fg <- factor(g)
  prior <- table(fg)
  prior <- prior/sum(prior)
  cnames <- levels(fg)
  g <- as.numeric(fg)
  J <- length(cnames)
  iswt <- FALSE
  if (missing(weights)) 
    dp <- table(g)/n
  else {
    weights <- (n * weights)/sum(weights)
    dp <- tapply(weights, g, sum)/n
    iswt <- TRUE
  }
  if(is.null(list.length) || is.null(genelist))
    stop("A list.length and genelist must be given for the top.genes option.\n", call. = FALSE)
  if (missing(theta)) 
    theta <- contr.treatment(J)
  theta <- contr.fda(dp, theta)
  Theta <- theta[g, , drop = FALSE]
  method <- match.arg(method)
  fit <- switch(method,
                pls = svdpls1c(scale(x, center = TRUE, scale = FALSE),
                  scale(Theta, center = TRUE, scale = FALSE), dimension),
                pcr = svdr(x, Theta, dimension + 1),
                ridge = gen.ridge(x, Theta, weights, lambda = 1, ...))
  genes <- vector("list", length = dim(theta)[2])
  for(i in seq(along = genes)){
    ord <- order(abs(coef(fit)[,i]), decreasing = TRUE)
    genes[[i]] <- genelist[ord][1:list.length]
  }
  gps <- vector("list", length(unique(fg)))
  newx <- matrix(NA, nc = dim(x)[2], nr = dim(x)[1])
  counts  <- vector("list", dim(theta)[2])
  for(k in seq(along = counts)){
    counts[[k]] <- matrix(NA, nc = B, nr = list.length)
    row.names(counts[[k]]) <- genes[[k]]
  }
  nam <- paste(levels(y)[as.numeric(dimnames(contr.treatment(J))[[2]])], "vs",
               levels(y)[1])
  names(counts) <- nam
  for(i in unique(fg)) gps[[i]] <- which(fg == unique(fg)[i])
  for(i in 1:B){
    for(j in seq(along = gps)){
      samp <- sample(gps[[j]], length(gps[[j]]), TRUE)
      newx[gps[[j]],] <- t(apply(x[samp,], 1, jitter))
    }
    fit.tmp <-  switch(method,
                       pls = svdpls1c(scale(newx, center = TRUE, scale = FALSE),
                         scale(Theta, center = TRUE, scale = FALSE), dimension),
                       pcr = svdr(newx, Theta, dimension + 1),
                       ridge = gen.ridge(newx, Theta, weights, lambda = 1, ...))
    for(k in seq(along = counts)){
      ord <- order(abs(coef(fit.tmp)[,k]), decreasing = TRUE)
      counts[[k]][,i] <- genes[[k]] %in% genelist[ord][1:list.length]
    }
  }
  counts <- lapply(counts, function(x) rowSums(x)/B)
  counts <- lapply(counts, as.data.frame)
  counts
}

pdmClass.cv <- function(Y, X, method = c("pls", "pcr", "ridge")){
  out <- vector()
  out <- factor(out, levels = levels(Y))
  for(i in seq(along = Y)){
    tmp <- pdmClass(Y[-i] ~ X[-i,], method = method)
    out[i] <- predict(tmp, X[i, , drop = FALSE])
  }
  out
}
    
