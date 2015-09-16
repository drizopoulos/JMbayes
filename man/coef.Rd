\name{coef}

\alias{coef.JMbayes}
\alias{fixef.JMbayes}
\alias{confint.JMbayes}

\title{Estimated Coefficients and Confidence Intervals for Joint Models}

\description{
  Extracts estimated coefficients and confidence intervals from fitted joint models.
}

\usage{
\method{coef}{JMbayes}(object, process = c("Longitudinal", "Event"), \dots)

\method{fixef}{JMbayes}(object, process = c("Longitudinal", "Event"), \dots)

\method{confint}{JMbayes}(object, parm = c("all", "Longitudinal", "Event"), \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMbayes}.}
  \item{process}{for which submodel (i.e., linear mixed model or survival model) to extract the estimated 
    coefficients.}
  \item{parm}{for which submodel (i.e., linear mixed model or survival model) to extract credible intervals.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  When \code{process = "Event"} both methods return the same output. However, for \code{process = "Longitudinal"},
  the \code{coef()} method returns the subject-specific coefficients, whereas \code{fixef()} only the fixed effects.
}

\value{A numeric vector or a matrix of the estimated parameters or confidence intervals for the fitted model.}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{ranef.JMbayes}}, \code{\link{jointModelBayes}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# joint model fit
fitJOINT <- jointModelBayes(fitLME, fitCOX, 
    timeVar = "obstime")

# fixed effects for the longitudinal process
fixef(fitJOINT)

# fixed effects + random effects estimates for the longitudinal 
# process
coef(fitJOINT)

# fixed effects for the event process
fixef(fitJOINT, process = "Event")
coef(fitJOINT, process = "Event")
}
}

\keyword{ methods }
