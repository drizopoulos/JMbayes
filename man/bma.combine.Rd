\name{bma.combine}

\alias{bma.combine}

\title{ Combines Predictions for Bayesian Model Averaging }

\description{
    Combines estimated survival probabilities or predictions for longitudinal responses.
}

\usage{
bma.combine(\dots, JMlis = NULL, weights = NULL)
} 

\arguments{
  \item{\dots}{objects inheriting from class \code{survfit.JMbayes} or \code{predict.JMbayes}.}
  \item{JMlis}{a list of \code{survfit.JMbayes} or \code{predict.JMbayes} objects.}
  \item{weights}{a numeric vector of weights to be applied in each object.}
}

\value{
  an object of class \code{survfit.JMbayes} or \code{predict.JMbayes}.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D., Hatfield, L., Carlin, B. and Takkenberg, J. (2014). Combining dynamic predictions from joint models for longitudinal and time-to-event data using Bayesian model averaging. \emph{Journal of the American Statistical Association} 
\bold{109}, 1385--1397.
}


\examples{
\dontrun{
# we construct the composite event indicator (transplantation or death)
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# we fit two joint models using splines for the subject-specific 
# longitudinal trajectories and a spline-approximated baseline
# risk function; the first one with the current value parameterization
# and the other with the shared random effects parameterization
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
              random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

jointFit1 <- jointModelBayes(lmeFit, survFit, timeVar = "year")
jointFit2 <- update(jointFit1, param = "shared-RE")

# we compute survival probabilities for Subject 2 with 
# different weights
ND <- pbc2[pbc2$id == 2, ] # the data of Subject 2
survPreds1 <- survfitJM(jointFit1, newdata = ND, weight = 0.4)
survPreds2 <- survfitJM(jointFit2, newdata = ND, weight = 0.6)

survPreds.bma <- bma.combine(survPreds1, survPreds2)
survPreds.bma
plot(survPreds.bma)
}
}

\keyword{methods}

