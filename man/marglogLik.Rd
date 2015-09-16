\name{marglogLik}

\alias{marglogLik}

\title{ Calculates Marginal Subject-specific Log-Likelihood Contributions }

\description{
    This function computes marginal subject-specific log-likelihood contributions based on a fitted joint model. The marginalization is done
    with respect to both the random effects and the parameters using a Laplace approximation.
}

\usage{
marglogLik(object, newdata, idVar = "id", method = "BFGS", control = NULL)
} 

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{method}{the \code{method} argument of \code{optim()}.}
  \item{control}{the \code{control} argument of \code{optim()}.}
}

\value{
  a numeric vector of marginal log-likelihood contributions.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
\dontrun{
# we construct the composite event indicator (transplantation or death)
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# we fit a joint model using splines for the subject-specific 
# longitudinal trajectories and a spline-approximated baseline
# risk function
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
    random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

# we compute marginal log-likelihood contribution for Subject 2
ND <- pbc2[pbc2$id == 2, ] # the data of Subject 2
marglogLik(jointFit, newdata = ND)
}
}

\keyword{methods}

