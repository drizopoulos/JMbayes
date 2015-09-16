\name{ranef}

\alias{ranef.JMbayes}

\title{Random Effects Estimates for Joint Models}

\description{
Extracts the random effects estimates from a fitted joint model.
}

\usage{
\method{ranef}{JMbayes}(object, postVar = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMbayes}.}
  \item{postVar}{logical; if \code{TRUE} the variance of the posterior distribution is also returned.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{
  a numeric matrix with rows denoting the individuals and columns the random effects (e.g., intercepts, slopes, etc.). 
  If \code{postVar = TRUE}, the numeric matrix has an extra attribute ``postVar".
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\seealso{\code{\link{coef.JMbayes}}, \code{\link{jointModelBayes}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# joint model fit, under the (default) Weibull model
fitJOINT <- jointModelBayes(fitLME, fitSURV, timeVar = "year")
ranef(fitJOINT)
}
}

\keyword{methods}
