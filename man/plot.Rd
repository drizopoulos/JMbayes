\name{plot}

\alias{plot.JMbayes}

\title{MCMC Diagnostics for Joint Models}

\description{
  Produces MCMC diagnostics plots.
}

\usage{
\method{plot}{JMbayes}(x, which = c("trace", "autocorr", "density", "CPO", "weightFun"), 
    param = c("betas", "sigma", "D", "gammas", "alphas", "Dalphas", 
        "shapes", "Bs.gammas", "tauBs"), ask = TRUE, max.t = NULL, 
    from = 0, \dots)
}

\arguments{
  \item{x}{an object inheriting from class \code{JMbayes}.}
  \item{which}{which types of plots to produce.}
  \item{param}{for which parameter to produce the MCMC diagnostic plots; default is for all parameters.}
  \item{ask}{logical, if \code{TRUE} the user is asked for input, before a new figure is drawn.}
  \item{max.t}{numeric scalar; up to which time point to plot the weight function, default is up to the
	third quantile of the observed event times.}
  \item{from}{numeric scalar; from which time point to start plotting the weight function, default is zero.}
  \item{\dots}{additional arguments; currently none is used.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\seealso{\code{\link{jointModelBayes}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)

# survival regression fit
fitSURV <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# joint model fit, under the (default) Weibull model
fitJOINT <- jointModelBayes(fitLME, fitSURV, timeVar = "year")

plot(fitJOINT)
}
}

\keyword{methods}
