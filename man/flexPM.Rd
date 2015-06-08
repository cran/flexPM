\name{flexPM}
\alias{flexPM}
\alias{cflexPM}
\alias{ctflexPM}
\title{Flexible Modeling of Survival Data}
\description{
Fit a flexible parametric regression model to possibly right-censored, left-truncated data.
}
\usage{
flexPM(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6)
cflexPM(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6)
ctflexPM(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6)
}
\arguments{
  \item{formula}{
  an object of class \dQuote{\code{\link{formula}}} (or one that can be coerced to that class): 
  a symbolic description of the model to be fitted.
  See \sQuote{Details}.

}
  \item{data}{
  an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) 
  containing the variables in the model. If missing, the variables are taken from 
  \code{environment(formula)}.
}
  \item{weights}{optional weights to be used during the fitting process.}
  \item{df}{the degrees of freedom of the B-spline that describes \eqn{u(T)} (see \sQuote{Details}).}
  \item{degree}{the degree of the polynomial of the B-spline that describes \eqn{u(T)}.}
  \item{knots}{the \emph{internal} knots of the B-spline that describes \eqn{u(T)}.}
  \item{maxit}{maximum number of iterations of the Newton-Raphson algorithm. If missing, 
	 a default value based on the number of free parameters is used.}
   \item{tol}{tolerance for the Newton-Raphson algorithm.}
}

\details{
This function fits a flexible parametric model to possibly right-censored, left-truncated outcomes,
usually survival data. Right censoring occurs when instead of some outcome \eqn{T}, one can only observe 
\eqn{Y = \min(T,C)}{Y = min(T,C)} and \eqn{d = I(T \le C)}, the indicator of failure. 
Left truncation occurs when \eqn{Y} is only observed subject to \eqn{Y \ge Z}.
Typically, \eqn{Z} is the time at enrollement of subjects with a disease. Those who died 
before enrollment (\eqn{Y < Z}) cannot be observed, thus generating selection bias.


The \code{formula} should be of the form 
\itemize{
\item{\code{y ~ x}, when using \command{flexPM} with non censored, non truncated data;}
\item{\code{cbind(y,d) ~ x}, when using \command{cflexPM} with censored data;}
\item{\code{cbind(y,z,d) ~ x}, when using \command{ctflexPM} with censored, truncated data.}
}
In the above, \code{x} is a set of predictors, \code{y} is the response variable,
\code{z} truncation times, and \code{d} the indicator of failure (1 = event, 0 = censored).
If \code{d} is omitted, it will be assumed to be always \kbd{1}.
Non-truncated observations are allowed in \command{ctflexPM} by letting \code{z = -Inf}.

The model is fitted by assuming that some transformation \eqn{u(T)} of \eqn{T} follows a logistic distribution (see \code{\link{Logistic}}),
conditional on a vector \eqn{x} of observed covariates.
This means that the following cumulative distribution function is used:
\deqn{F_T(t \mid x) = \frac{1}{1 + e^{-\frac{u(t) - m(x)}{s(x)}}}}{F(t|x) = 1/exp\{-[u(t) - m(x)]/s(x)\}}
where \eqn{m(x)} and \eqn{\log s(x)}{log s(x)} are modeled as specified by
\code{formula}, while \eqn{u(T)} is defined as a B-spline function built via \command{spline.des} (see \code{\link{bs}}).
You can choose the degrees of freedom \kbd{df} and the \kbd{degree} of the spline basis. Estimation is carried out by Newton-Raphson algorithm. The model parameters are (a) the coefficients
describing the effect of covariates \eqn{x} on \eqn{m(x)}
and \eqn{\log s(x)}{log s(x)}, 
and (b) the coefficients of the B-spline basis that defines the unknown transformation \eqn{u(\cdot)}{u(.)}.
}
\value{
An object of class \dQuote{\code{flexPM}}, which is a list with the following items:
 \item{converged}{logical value indicating the convergence status of the algorithm.}
 \item{n.it}{the number of iterations.}
 \item{n}{the number of observations.}
 \item{n.free.par}{the number of free parameters in the model.}
  \item{loglik}{the values of the log-likelihood at convergence.}
  \item{AIC, BIC}{the Akaike and Bayesian information criterion.}
\item{mf}{the used model frame.}
\item{call}{the matched call.}
The model parameters are returned as attributes of \code{mf} and are not
user-level objects. The fit is only intended for prediction:
use \code{\link{predict.flexPM}}.
}
\note{
The choice of the logistic distribution is only due to computational convenience and
any other family of distributions could be used. Note, however, that distributional 
assumptions are made on \eqn{u(T)}, which is itself \emph{estimated}, 
and not directly on \eqn{T}: this permits achieving a great flexibility 
and enables fitting asymmetric or multi-modal distributions. Which distribution
is used as \dQuote{kernel} is quite irrelevant.
This flexible parametric approach generally outoperforms fully
nonparametric estimators like the local Kaplan-Meier, at a cost of a relatively small bias. 


The model is intended to be used for prediction and not for inference. The hessian matrix is not returned.
The number of free parameters is \code{df + 2*ncol(x) - 1}, and not \code{df + 2*ncol(x)}, because
the scale of \eqn{u(\cdot)}{u(.)} and that of the logistic distribution are exchangeable and thus one 
coefficient of \eqn{u(\cdot)}{u(.)} is constrained to be \kbd{1}. 
Note that with \kbd{df = 1} and \kbd{degree = 1}, \eqn{u(\cdot)}{u(.)} is linear,
i.e., you are actually fitting a logistic distribution.

Finally, be aware that you cannot
remove the intercept from the model: \code{~ -1} in  \code{formula} is ignored.
}

\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\seealso{
\code{\link{predict.flexPM}}
}
\examples{
# Simulated data from a normal distribution 

set.seed(1234); n <- 1000
x1 <- rnorm(n)
x2 <- runif(n)


# non-censored, non-truncated data

t <- rnorm(n, 2 + 3*x1, 1 + x2) # time variable
m1 <- flexPM(t ~ x1 + x2)


# right-censored data

c <- rnorm(n,3,3) # censoring variable
y <- pmin(t,c)    # observed outcome
d <- (t <= c)     # 1 = observed, 0 = censored
m2 <- cflexPM(cbind(y,d) ~ x1 + x2)


# right-censored, left-truncated data

z <- rnorm(n,-3,3) # truncating variable
w <- which(y > z)  # only observe if y > z
y <- y[w]; d <- d[w]; z <- z[w]; x1 <- x1[w]; x2 <- x2[w]
m3 <- ctflexPM(cbind(y,z,d) ~ x1 + x2)


# m1, m2, m3 are not intended to be interpreted.
# Use the loglik, AIC and BIC for model selection 
# (choice of df, inclusion/exclusion of covariates).
# Use predict() to obtain predictions.

# Note that the following are identical:
# flexPM(y ~ x1 + x2)
# cflexPM(cbind(y,1) ~ x1 + x2)
# ctflexPM(cbind(y,-Inf,1) ~ x1 + x2)
}


\keyword{models}
\keyword{survival}



