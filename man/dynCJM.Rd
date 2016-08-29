\name{dynCJM}

\alias{dynCJM}
\alias{dynCJM.JMbayes}

\title{A Dynamic Discrimination Index for Joint Models}

\description{ 
  This function computes a dynamic discrimination index for joint models based on a weighted average of time-dependent AUCs.
}

\usage{
dynCJM(object, newdata, Dt, \dots)

\method{dynCJM}{JMbayes}(object, newdata, Dt, idVar = "id", t.max = NULL, 
    simulate = FALSE, M = 100, weightFun = NULL, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{Dt}{a numeric scalar denoting the time frame within which the occurrence of events is of interest.}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{t.max}{a numeric scalar denoting the time maximum follow-up time up to which the dynamic discrimination index is to be calculated.
    If \code{NULL}, it is set equal to \code{max(Time) + 1e-05} where \code{Time} denotes the observed event times.}
  \item{simulate}{logical; if \code{TRUE}, a Monte Carlo approach is used to estimate survival probabilities. If \code{FALSE},
    a first order estimator is used instead. See \code{\link{survfitJM}} for mote details.}
  \item{M}{a numeric scalar denoting the number of Monte Carlo samples; see \code{\link{survfitJM}} for mote details.}
  \item{weightFun}{a function of two arguments the first denoting time and the second the length of the time frame of interest, i.e., \code{Dt}.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
(\bold{Note:} The following contain some math formulas, which are better viewed in the pdf version 
of the manual accessible at \url{http://cran.r-project.org/package=JMbayes}.)

Function \code{dynC} computes the following discrimination index
\deqn{\mbox{C}_{dyn}^{\Delta t} = \int_0^{t_{max}} \mbox{AUC}(t, \Delta t) \, 
\mbox{Pr} \{ {\cal E}(t, \Delta t) \} \; dt \Big / \int_0^{t_{max}} \mbox{Pr} \{ {\cal E}(t, \Delta t) \} \; dt,}
where \deqn{\mbox{AUC}(t, \Delta t) = \mbox{Pr} \bigl [ \pi_i(t + \Delta t \mid t) < 
\pi_j(t + \Delta t \mid t) \mid \{ T_i^* \in (t, t + \Delta t] \} \cap \{ T_j^* > t + \Delta t \} \bigr ],}
and \deqn{{\cal E}(t, \Delta t) = \bigl [ \{ T_i^* \in (t, t + \Delta t] \} \cap \{ T_j^* > t + 
\Delta t \} \bigr ],} with \eqn{i} and \eqn{j} denote a randomly selected pair subjects, and
\eqn{\pi_i(t + \Delta t \mid t)} and \eqn{\pi_j(t + \Delta t \mid t)} denote the conditional survival probabilities calculated by
\code{\link{survfitJM}} for these two subjects, for different time windows \eqn{\Delta t} specified by argument \code{Dt}. 
The upper limit of integral in specified by argument \code{t.max}. The integrals in the numerator and denominator 
are approximated using a 15-point Gauss-Kronrod quadrature rule.
}

\value{
  A list of class \code{dynCJM} with components:
  \item{dynC}{a numeric scalar denoting the dynamic discrimination index.}
  \item{times}{a numeric vector of time points at which the AUC was calculated.}
  \item{AUCs}{a numeric vector of the estimated AUCs at the aforementioned time points.}
  \item{weights}{a numeric vector of the estimated weights at the aforementioned time points.}
  \item{t.max}{a copy of the \code{t.max} argument.}
  \item{Dt}{a copy of the \code{Dt} argument.}
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


\seealso{\code{\link{survfitJM}}, \code{\link{aucJM}}, \code{\link{jointModelBayes}}}

\examples{
\dontrun{
# we construct the composite event indicator (transplantation or death)
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# we fit the joint model using splines for the subject-specific 
# longitudinal trajectories and a spline-approximated baseline
# risk function
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
              random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

# dynamic discrimination index up to year 10 using a two-year interval 
dynCJM(jointFit, pbc2, Dt = 2, t.max = 10)
}
}

\keyword{methods}
