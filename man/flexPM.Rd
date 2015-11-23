\name{flexPM}
\alias{flexPM}
\title{Flexible Modeling of Survival Data}
\description{
Fit a flexible parametric regression model to possibly right-censored, left-truncated data.
}
\usage{
flexPM(formula, data, weights, df = 3, degree = 3, knots, maxit, tol = 1e-6)
}
\arguments{
  \item{formula}{
  an object of class \dQuote{\code{\link{formula}}}: 
  a symbolic description of the model to be fitted. The response must be a survival object 
  as returned by the \code{Surv} function.
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
Left truncation occurs when \eqn{Y} is only observed subject to \eqn{Y > Z}.
Typically, \eqn{Z} is the time at enrollement of subjects with a disease. Those who died 
before enrollment (\eqn{Y < Z}) cannot be observed, thus generating selection bias.


The \code{formula} must be of the form 
\itemize{
\item{\code{Surv(y,d) ~ x}, with censored data;}
\item{\code{Surv(z,y,d) ~ x}, with censored, truncated data;}
\item{\code{Surv(y) ~ x}, is also allowed and denotes non-censored, non-truncated data.}
}
In the above, \code{x} is a set of predictors, \code{y} is the response variable,
\code{z} truncation times (\code{z < y}), and \code{d} the indicator of failure (1 = event, 0 = censored).

In \code{flexPM}, model fitting is implemented as follows.
First, the response variable is pre-trasformed using a smoothed version of
\kbd{y = qlogis(rank(y)/(n + 1))}. Second, parameter estimation is carried out on the transformed variable. 
Maximum likelihood estimators are computed via Newton-Raphson algorithm, using the following flexible distribution:
\deqn{F_T(t \mid x) = \frac{1}{1 + e^{-\frac{u(t) - m(x)}{s(x)}}}.}{F(t|x) = 1/(1 + exp\{-[u(t) - m(x)]/s(x)\}).}
In the above, \eqn{m(x)} and \eqn{\log s(x)}{log s(x)} are modeled as specified by
\code{formula}, while \eqn{u(\cdot)}{u(.)} is a B-spline function built via \command{spline.des} (see \code{\link{bs}}).
You can choose the degrees of freedom \kbd{df} and the \kbd{degree} of the spline basis. The model parameters are (a) the coefficients
describing the effect of covariates \eqn{x} on \eqn{m(x)}
and \eqn{\log s(x)}{log s(x)}, 
and (b) the coefficients of the B-spline basis that defines the unknown transformation \eqn{u(\cdot)}{u(.)}, on which suitable constraints are imposed to ensure monotonicity.
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
user-level objects. The model is intended to be used for prediction and not for inference. The hessian matrix 
is not returned. The number of free parameters is \code{df + 2*ncol(x) - 1}, 
and not \code{df + 2*ncol(x)}, because the scale of \eqn{u(\cdot)}{u(.)} and that of \eqn{F_T(t \mid x)}{F(t|x)} 
are exchangeable and thus one coefficient of \eqn{u(\cdot)}{u(.)} is constrained to be \kbd{1}. 

The accessor functions \code{summary}, \code{nobs}, \code{logLik}, \code{AIC}, and \code{BIC} can be used 
to extract information from the model. The fit is only intended for prediction: use \code{\link{predict.flexPM}}.
}
\note{
The model is fitted assuming that an unknown transformation of the response 
variable follows a Logistic distribution. The choice of the \dQuote{kernel} distribution is only due to computational convenience and does not reflect any prior belief. Provided that \eqn{u(\cdot)}{u(.)} is sufficiently flexible, asymmetric or multi-modal distributions can be fitted. Pre-transforming the response variable (added in \pkg{flexPM 2.0}) removes the outliers, generates a more symmetric distribution, and frequently permits achieving a good fit using fewer knots for \eqn{u(\cdot)}{u(.)}.

This flexible parametric approach generally outoperforms fully nonparametric estimators like local Kaplan-Meier, 
at a cost of a relatively small bias. 
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\seealso{
\code{\link{predict.flexPM}}
}
\examples{
# Simulated data from a normal distribution 

n <- 1000
x1 <- rnorm(n)
x2 <- runif(n)


# non-censored, non-truncated data

t <- rnorm(n, 2 + 3*x1, 1 + x2) # time variable
m1 <- flexPM(Surv(t) ~ x1 + x2)


# right-censored data

c <- rnorm(n,3,3) # censoring variable
y <- pmin(t,c)    # observed outcome
d <- (t <= c)     # 1 = observed, 0 = censored
m2 <- flexPM(Surv(y,d) ~ x1 + x2)


# right-censored, left-truncated data

z <- rnorm(n,-3,3) # truncating variable
w <- which(y > z)  # only observe if y > z
y <- y[w]; d <- d[w]; z <- z[w]; x1 <- x1[w]; x2 <- x2[w]
m3 <- flexPM(Surv(z,y,d) ~ x1 + x2)

################################################################

# m1, m2, m3 are not intended to be interpreted.
# Use predict() to obtain predictions.

# Note that the following are identical:
# summary(flexPM(Surv(y) ~ x1 + x2))
# summary(flexPM(Surv(y, rep(1,length(y))) ~ x1 + x2))
# summary(flexPM(Surv(rep(-Inf,length(y)), y, rep(1,length(y))) ~ x1 + x2))

################################################################

# Use the logLik, AIC and BIC for model selection 
# (choice of df, inclusion/exclusion of covariates)

models <- list(
  flexPM(Surv(z,y,d) ~ x1 + x2, df = 1, degree = 1),
  flexPM(Surv(z,y,d) ~ x1 + x2, df = 3),
  flexPM(Surv(z,y,d) ~ x1 + x2 + I(x1^2) + I(x2^2), df = 1, degree = 1),
  flexPM(Surv(z,y,d) ~ x1 + x2 + I(x1^2) + I(x2^2), df = 3),
  flexPM(Surv(z,y,d) ~ x1 * x2, df = 5)
)

my_final_model <- models[[which.min(sapply(models, function(x) x$AIC))]]
summary(my_final_model)

}


\keyword{models}
\keyword{survival}



