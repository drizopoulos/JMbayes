\name{xtable}

\alias{xtable.JMbayes}

\title{xtable Method from Joint Models.}

\description{
  produces a LaTeX table with the results of a joint model using package xtable.
}

\usage{
\method{xtable}{JMbayes}(x, caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL, which = c("all", "Longitudinal", "Event"),
    varNames.Long = NULL, varNames.Event = NULL, \dots)
}

\arguments{
  \item{x}{an object inheriting from class \code{JMbayes}.}
  \item{caption}{the \code{caption} argument of \code{xtable()}.}
  \item{label}{the \code{label} argument of \code{xtable()}.}
  \item{align}{the \code{align} argument of \code{xtable()}.}
  \item{digits}{the \code{digits} argument of \code{xtable()}.}
  \item{display}{the \code{display} argument of \code{xtable()}.}
  \item{which}{a character string indicating which results to include in the
      LaTeX table. Options are all results, the results of longitudinal
      submodel or the results of the survival submodel.}
  \item{varNames.Long}{a character vector of the variable names for the
      longitudinal submodel.}
  \item{varNames.Event}{a character vector of the variable names for the
      survival submodel.}
  \item{\dots}{additional arguments; currently none is used.}
}

\value{A LaTeX code chunk with the results of the joint modeling analysis.}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{jointModelBayes}}}

\examples{
\dontrun{
prothro$t0 <- as.numeric(prothro$time == 0)
lmeFit <- lme(pro ~ treat * (time + t0), random = ~ time | id, data = prothro)
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)
jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "time")
    
if (require("xtable")) {
    xtable:::xtable(jointFit, math.style.negative = TRUE)
}
}
}

\keyword{ methods }
