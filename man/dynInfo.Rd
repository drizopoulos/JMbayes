\name{dynInfo}

\alias{dynInfo}

\title{Dynamic Information of an Extra Longitudinal Measurement}

\description{ 
  Using the available longitudinal information up to a particular time point, this 
  function computes an estimate of the information we again by obtaining an extra
  longitudinal measurement.
}

\usage{
dynInfo(object, newdata, Dt, K = 5, M = 500, idVar = "id", 
    simulateFun = function (eta, scale) rnorm(length(eta), eta, scale), 
    seed = 1L)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for
  the subject for whom we wish to plan the next measurement. The names of the variables in
  this data frame must be the same as in the data frames that were used to fit the linear 
  mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()}
  that were supplied as the two first argument of \code{\link{jointModelBayes}}. In 
  addition, this data frame should contain a variable that identifies the subject 
  (see also argument \code{idVar}).}
  \item{Dt}{numeric scalar denoting the length of the time interval to search for the 
  optimal time point of the next measurement, i.e., the interval is \eqn{(t, t + Delta t]}
  with \eqn{Delta t} given by \code{Dt}.}
  \item{K}{numeric scalar denoting the number of time points to cosider in the interval 
  \eqn{(t, t + Delta t]}.}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the subject.}
  \item{simulateFun}{a function based on which longitudinal measurement can be simulated. 
  This should have as a main argument the variable \code{eta} that denotes the 
  subject-specific linear predictor from the mixed model, and possibly a \code{scale} 
  parameter.}
  \item{M}{a numeric scalar denoting the number of Monte Carlo samples.}
  \item{seed}{a numeric scalar}
}

\details{
This functions computes the following posterior predictive distribution
\deqn{E_{Y} [ E_{T^* | Y} (\log p (T_j^* \mid T_j^* > u, \{ Y_j(t), y_j(u) \}, 
 D_n \bigr )) ],} where \eqn{T_j^*} denotes the time-to-event for subject \eqn{j} for
 whom we wish to plan the next visit, \eqn{Y_j(t)} the available longitudinal measurements
 of this subject up to time \eqn{t}, \eqn{y_j(u)} the future longitudinal measurement we
 wish to plan at time \eqn{u > t}, and \eqn{D_n} the data set that was used to fit the 
 joint model.
}

\value{
  A list with components:
  \item{summary}{a numeric matrix with first column the time points at which the
  longitudinal measurement is hypothetically taken, second column the estimated 
  information we gain by obtaining the measurement, and third column the estimated 
  cumulative risk of an event up to the particular time point denoted in the first column.}
  \item{full.results}{a numeric matrix with columns representing the time points, rows
  the Monte Carlo samples, and entries the value of log posterio predictive density.}
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{survfitJM}}, \code{\link{jointModelBayes}}}

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

dynInfo(jointFit, newdata = pbc2[pbc2$id == 2, ], Dt = 5)[[1]]
}
}

\keyword{methods}
