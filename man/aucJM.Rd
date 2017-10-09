\name{aucJM}

\alias{aucJM}
\alias{rocJM}
\alias{aucJM.JMbayes}
\alias{aucJM.mvJMbayes}
\alias{rocJM.JMbayes}
\alias{rocJM.mvJMbayes}


\title{Time-Dependent ROCs and AUCs for Joint Models}

\description{ 
  Using the available longitudinal information up to a starting time point, this function computes an estimate 
  of the ROC and the AUC at a horizon time point based on joint models.
}

\usage{
aucJM(object, newdata, Tstart, \dots)

\method{aucJM}{JMbayes}(object, newdata, Tstart, Thoriz = NULL, 
    Dt = NULL, idVar = "id", simulate = FALSE, M = 100, \dots)

\method{aucJM}{mvJMbayes}(object, newdata, Tstart, Thoriz = NULL, 
    Dt = NULL, idVar = "id", M = 100, \dots)

rocJM(object, newdata, Tstart, \dots)

\method{rocJM}{JMbayes}(object, newdata, Tstart, Thoriz = NULL, 
    Dt = NULL, idVar = "id", simulate = FALSE, M = 100, \dots)

\method{rocJM}{mvJMbayes}(object, newdata, Tstart, Thoriz = NULL, 
    Dt = NULL, idVar = "id", M = 100, \dots)
    
predict_eventTime(object, newdata, cutpoints, \dots)

\method{predict_eventTime}{mvJMbayes}(object, newdata, cut_points, 
    idVar = "id", M = 500L, low_percentile = 0.025, \ldots)

find_thresholds(object, newdata, Dt, \dots)

\method{find_thresholds}{mvJMbayes}(object, newdata, Dt, idVar = "id", 
    M = 200L, n_cores =  max(1, parallel::detectCores() - 2), \ldots)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMbayes} or \code{mvJMbayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; 
    \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL}
    is set equal to \code{Tstart + Dt}.}
  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{simulate}{logical; if \code{TRUE}, a Monte Carlo approach is used to estimate survival probabilities. If \code{FALSE},
    a first order estimator is used instead. See \code{\link{survfitJM}} for mote details.}
  \item{M}{a numeric scalar denoting the number of Monte Carlo samples; see \code{\link{survfitJM}} for mote details.}
  \item{cut_points}{a numeric matrix with first column time-points followed by other columns of 
  optimal cut-points from an ROC curve.}
  \item{low_percentile}{a numeric value indicating the percentile based on which it will
  be judged whether the spread of the posterior predictive distribution is too large.}
  \item{n_cores}{an integer indicating the number of cores to use for parallel computing.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
Based on a fitted joint model (represented by \code{object}) and using the data supplied in argument \code{newdata}, this function
computes the following estimate of the AUC: \deqn{\mbox{AUC}(t, \Delta t) = \mbox{Pr} \bigl [ \pi_i(t + \Delta t \mid t) < 
\pi_j(t + \Delta t \mid t) \mid \{ T_i^* \in (t, t + \Delta t] \} \cap \{ T_j^* > t + \Delta t \} \bigr ],}
with \eqn{i} and \eqn{j} denote a randomly selected pair of subjects, and
\eqn{\pi_i(t + \Delta t \mid t)} and \eqn{\pi_j(t + \Delta t \mid t)} denote the conditional survival probabilities calculated by
\code{\link{survfitJM}} for these two subjects, for different time windows \eqn{\Delta t} specified by argument \code{Dt} using
the longitudinal information recorded up to time \code{t =} \code{Tstart}.

The estimate of \eqn{\mbox{AUC}(t, \Delta t)} provided by \code{aucJM()} is in the spirit of Harrell's 
\eqn{c}-index, that is for the comparable subjects (i.e., the ones whose observed event times can be ordered), we 
compare their dynamic survival probabilities calculated by \code{\link{survfitJM}}. For the subjects who due to 
censoring we do not know if they are comparable, they contribute in the AUC with the probability that they would
have been comparable.
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

\references{
Antolini, L., Boracchi, P., and Biganzoli, E. (2005). A time-dependent discrimination index
for survival data. \emph{Statistics in Medicine} \bold{24}, 3927--3944.

Harrell, F., Kerry, L. and Mark, D. (1996). Multivariable prognostic models: issues in 
developing models, evaluating assumptions and adequacy, and measuring and reducing errors.
\emph{Statistics in Medicine} \bold{15}, 361--387.

Heagerty, P. and Zheng, Y. (2005). Survival model predictive accuracy and ROC curves. 
\emph{Biometrics} \bold{61}, 92--105.

Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for 
longitudinal and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.
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

# AUC using data up to year 5 with horizon at year 8 
aucJM(jointFit, pbc2, Tstart = 5, Thoriz = 8)

plot(rocJM(jointFit, pbc2, Tstart = 5, Thoriz = 8))
}
}

\keyword{methods}
