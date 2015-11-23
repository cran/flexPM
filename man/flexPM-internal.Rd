\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
BS(x, df, knots, degree, Boundary.knots = range(x), deriv = 0)
flexPM.newton(start, f, tol = 1e-05, maxit = 200, \ldots)
L.estimator(z,y,d,x,weights, df,degree,knots, theta, maxit = 100, tol = 1e-6, type)
L.loglik(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
L.wloglik(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
L.loglik.cens(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
L.wloglik.cens(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
L.loglik.cens.trunc(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
L.wloglik.cens.trunc(theta, Bz, By, B1y, x, d, d1, dd, weights, notrunc, deriv = 2)
\method{print}{flexPM}(x, \ldots)
\method{summary}{flexPM}(object, \ldots)
\method{print.summary}{flexPM}(x, \ldots)
\method{nobs}{flexPM}(object, \ldots)
\method{logLik}{flexPM}(object, \ldots)
}
\keyword{internal}

