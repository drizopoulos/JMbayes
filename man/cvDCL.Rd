\name{cvDCL}

\alias{cvDCL}

\title{Dynamic Information}

\description{ 
  Using the available longitudinal information up to a starting time point, this function computes an estimate 
  of the cross-entropy function based on joint models.
}

\usage{
cvDCL(object, newdata, Tstart, idVar = "id", M = 300L, seed = 123L)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{Tstart}{a numeric scalar indicating at which time to compute the cross-entropy.}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
   \item{M}{a numeric scalar denoting the number of Monte Carlo samples.}
  \item{seed}{a numeric scalar}
}

\details{
This function calculates an estimate of the cross-entropy at any time point \eqn{t} (given in \code{Tstart}) that can be
used to identify the joint model that best predicts future event giver survival up to \eqn{t}.
}

\value{
  A list of class \code{aucJM} with components:
  \item{auc}{a numeric scalar denoting the estimated prediction error.}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}  
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{survfitJM}}, \code{\link{dynCJM}}, \code{\link{jointModelBayes}}}

\examples{
\dontrun{
# we construct the composite event indicator (transplantation or death)
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# we fit the joint model using splines for the subject-specific 
# longitudinal trajectories and a spline-approximated baseline
# risk function
lmeFit <- lme(log(serBilir) ~ ns(year, 3),
    random = list(id = pdDiag(form = ~ ns(year, 3))), data = pbc2)
survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

cvDCL(jointFit, Tstart = 5)
}
}

\keyword{methods}
