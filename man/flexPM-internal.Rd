\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
BS(x, df, knots, degree, Boundary.knots = range(x), deriv = 0)
flexPM.internal(mf,cl, df, degree, knots, maxit, tol, type)
flexPM.newton(start, f, tol = 1e-05, maxit = 200, \ldots)
L.estimator(y,d,x,weights, df,degree,knots, theta, z, maxit = 100, tol = 1e-6, type)
L.loglik(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
L.wloglik(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
L.loglik.cens(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
L.wloglik.cens(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
L.loglik.cens.trunc(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
L.wloglik.cens.trunc(theta, By, B1y, x, d, d1, dd, weights, Bz, notrunc, deriv = 2)
\method{print}{flexPM}(x, \ldots)
}
\keyword{internal}

