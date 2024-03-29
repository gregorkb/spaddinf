% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spaddinf.R
\name{spadd.presmth.Legr.cv}
\alias{spadd.presmth.Legr.cv}
\title{Choose tuning parameters for fitting the desparsified lasso presmoothing estimator with Legendre polynomials}
\usage{
spadd.presmth.Legr.cv(
  X,
  Y,
  d.pre,
  K = 1,
  n.lambda,
  n.eta,
  n.folds,
  plot = FALSE
)
}
\arguments{
\item{X}{the design matrix}

\item{Y}{the response vector}

\item{d.pre}{the number of intervals in which to divide the support of each covariate}

\item{K}{the order of the Legendre polynomials. E.g. \code{K=0} fits piecwise constant, \code{K=1} fits piecewise linear functions.}

\item{n.lambda}{the number of candidate lambda values}

\item{n.eta}{the number of candidate eta values}

\item{n.folds}{the number of crossvalidation folds}

\item{plota}{logical indicating whether crossvalidation output should be plotted}
}
\value{
a list with the chosen values of the tuning parameters
}
\description{
Choose tuning parameters for fitting the desparsified lasso presmoothing estimator with Legendre polynomials
}
\examples{
data <- data_gen(n = 200,q = 50,r = .9)

spadd.presmth.Legr.cv.out <- spadd.presmth.Legr.cv(X = data$X,
                                                   Y = data$Y,
                                                   d.pre = 10,
                                                   n.lambda = 25,
                                                   n.eta = 25,
                                                   n.folds = 5,
                                                   plot = TRUE)

spadd.presmth.Legr.out <- spadd.presmth.Legr(X = data$X,
                                             Y = data$Y,
                                             d.pre = 10,
                                             lambda = spadd.presmth.Legr.cv.out$cv.lambda,
                                             eta = spadd.presmth.Legr.cv.out$cv.eta,
                                             n.foi = 6)

plot_presmth_Legr(x = spadd.presmth.Legr.out,
                  true.functions = list( f = data$f,
                                         X = data$X))
}
