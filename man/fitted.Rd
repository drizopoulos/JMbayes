\name{fitted & residuals}

\alias{fitted.JMbayes}
\alias{residuals.JMbayes}


\title{Fitted Values and Residuals for Joint Models}

\description{
  Calculates fitted values for joint models.
}

\usage{
\method{fitted}{JMbayes}(object, 
    process = c("Longitudinal", "longitudinal", "Event", "event"), 
    type = c("Marginal", "marginal", "Subject", "subject"), nullY = FALSE, \dots)

\method{residuals}{JMbayes}(object, 
    process = c("Longitudinal", "longitudinal", "Event", "event"), 
    type = c("Marginal", "marginal", "Subject", "subject", 
             "Martingale", "martingale", "nullMartingale", "nullmartingale"),
    standardized = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{process}{for which model (i.e., linear mixed model or survival model) to calculate fitted values or residuals.}
  \item{type}{what type of fitted values or residuals to calculate. See \bold{Details}.}
  \item{nullY}{logical; if \code{TRUE} the association parameters that connect the longitudinal and event time process 
    are set to zero.}
  \item{standardized}{logical; if \code{TRUE} standardized residuals are calculated.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  For \code{process = "Longitudinal"}, let \eqn{X} denote the design matrix for the fixed effects \eqn{\beta}, and 
  \eqn{Z} the design matrix for the random effects \eqn{b}. Then for \code{type = "Marginal"} the fitted values are 
  \eqn{X \hat{\beta},} whereas for \code{type = "Subject"} they are \eqn{X \hat{\beta} + Z \hat{b}}, where \eqn{\hat{\beta}}
  and \eqn{\hat{b}} denote the corresponding posterior means for the fixed and random effects. The corresponding residuals 
  are calculated by subtracting the fitted values from the observed data \eqn{y}. If \code{type = "Subject"} and 
  \code{standardized = TRUE}, the residuals are divided by the estimated residual standard error.
  
  For \code{process = "Event"} function \code{fitted()} calculates the cumulative hazard function at each time point a longitudinal
  measurement has been recorded. If \code{nullY = TRUE}, then the cumulative hazard is calculated without the contribution of the
  longitudinal process. Function \code{residuals()} calculates the martingales residuals or the martingale residuals without the 
  contribution of the longitudinal process when \code{type = "nullMartingale"}.
}

\value{
  a numeric vector of fitted values or residuals.
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
\dontrun{
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
              random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)
jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

fitted(jointFit, process = "Event")
residuals(jointFit, type = "Subject", standardized = TRUE)
}
}

\keyword{methods}
