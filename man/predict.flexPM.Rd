\name{predict.flexPM}
\alias{predict.flexPM}
\title{
Prediction from Fitted Flexible Parametric Models}
\description{
Predicts the distribution function and simulates new data from a fitted model.
}
\usage{
\method{predict}{flexPM}(object, type = c("CDF", "QF", "sim"), newdata, p, \ldots)
}
\arguments{
  \item{object}{
an object of class \dQuote{\code{flexPM}} obtained with \code{\link{flexPM}}, 
\code{\link{cflexPM}}, or \code{\link{ctflexPM}}.
}
  \item{type}{
the type of prediction (see \sQuote{Details}).
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, 
the model frame of the object is used.
}
  \item{p}{
the order(s) of the quantile to be computed (for \code{type = "QF"})
}
\item{\ldots}{for future methods.}
}
\details{
\itemize{
\item{If \kbd{type = "CDF"} (the default), the fitted cumulative distribution function (\acronym{CDF})
and the corresponding probability density function (\acronym{PDF}) and survival
function (\acronym{SF}) are returned.}
\item{If \kbd{type = "QF"}, conditional quantiles of the specified order(s) are computed.}
\item{If \kbd{type = "sim"}, data are simulated from the fitted model.}
}
New data can be supplied: observe that for \kbd{type = "CDF"}, \code{newdata} must include
the values of the response variable, and not just the covariates.
}
\value{
\itemize{
  \item{If \kbd{type = "CDF"}, a named data frame with
  variables \code{log.f} (the fitted log-\acronym{PDF}), 
  \code{log.F} (the log-\acronym{CDF}) and \code{log.S} 
  (the log-\acronym{SF}).
}
\item{If \kbd{type = "QF"}, a named data frame containing the fitted conditional quantiles
  of the specified order(s) in different columns.
}
\item{
  If \kbd{type = "sim"}, a vector of simulated data from the fitted model. 
}
}
  All types of prediction are computed at \code{newdata}, if supplied, or at the observed data, otherwise.
}

\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\seealso{
\code{\link{flexPM}}
}

\examples{

# Using simulated data

set.seed(1111); n <- 1000
x <- runif(n)
t <- rnorm(n, 1 + x, 1 + x)
model <- flexPM(t ~ x + I(x^2)) 
# using polynomials (e.g. x^2) to achieve flexibility



# Prediction of the conditional cumulative distribution function (CDF)
# and the probability density function (PDF)

pred <- predict(model, type = "CDF") # predict the CDF and PDF

plot(pnorm(t, 1 + x, 1 + x), exp(pred$log.F))
abline(0,1, col = "green", lwd = 3) # true vs fitted CDF

plot(dnorm(t, 1 + x, 1 + x), exp(pred$log.f))
abline(0,1, col = "green", lwd = 3) # true vs fitted PDF



# Prediction of quantiles

predMe <- predict(model, type = "QF", p = 0.5) # predict the median
plot(x,t)
points(x, predMe$p0.5, col = "green") # fitted median
abline(1,1, col = "red", lwd = 3) # true median = 1 + x



# Simulate data from the fitted model

t.sim <- predict(model, type = "sim")
plot(sort(t.sim), sort(t)); abline(0,1)
# if the model is good, t and t.sim should have a similar distribution



######### Using new data #############################

newdata <- data.frame(t = c(0,1,2), x = c(0.1,0.5,0.9))
# note that new 't' is only needed for type = "CDF"

predict(model, type = "CDF", newdata = newdata)
predict(model, type = "QF", newdata = newdata, p = c(0.25,0.5,0.75))
predict(model, type = "sim", newdata = newdata)

}
\keyword{methods}
\keyword{distribution}

