\name{survfitJM}

\alias{survfitJM}
\alias{survfitJM.JMbayes}
\alias{survfitJM.mvJMbayes}

\title{Prediction in Joint Models}

\description{ 
  This function computes the conditional probability of surviving later times than the last observed time for which a
  longitudinal measurement was available.
}

\usage{
survfitJM(object, newdata, \dots)

\method{survfitJM}{JMbayes}(object, newdata, 
   type = c("SurvProb", "Density"), idVar = "id", 
   simulate = TRUE, survTimes = NULL, last.time = NULL, 
   LeftTrunc_var = NULL, M = 200L, 
   CI.levels = c(0.025, 0.975), log = FALSE, scale = 1.6, 
   weight = rep(1, nrow(newdata)), 
   init.b = NULL, seed = 1L, \dots)

\method{survfitJM}{mvJMbayes}(object, newdata, 
    survTimes = NULL, idVar = "id", last.time = NULL, 
    M = 200L, scale = 1.6, log = FALSE, 
    CI.levels = c(0.025, 0.975), seed = 1L, \dots)
}


\arguments{
  \item{object}{an object inheriting from class \code{JMBayes} or \code{mvJMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{type}{character string indicating what to compute, i.e., survival probabilities or the log conditional density.}
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{simulate}{logical; if \code{TRUE}, a Monte Carlo approach is used to estimate survival probabilities. If \code{FALSE},
    a first order estimator is used instead. (see \bold{Details})}
  \item{survTimes}{a numeric vector of times for which prediction survival probabilities are to be computed.}
  \item{last.time}{a numeric vector or character string. This specifies the known time at which each of the subjects in \code{newdata}
    was known to be alive. If \code{NULL}, then this is automatically taken as the last time each subject provided a longitudinal
    measurement. If a numeric vector, then it is assumed to contain this last time point for each subject. If a character string, then
    it should be a variable in the data frame \code{newdata}.}
  \item{LeftTrunc_var}{character string indicating the name of the variable in \code{newdata} that denotes the left-truncation
   time.}
  \item{M}{integer denoting how many Monte Carlo samples to use -- see \bold{Details}.}
  \item{CI.levels}{a numeric vector of length two that specifies which quantiles to use for the calculation of confidence interval for the 
    predicted probabilities -- see \bold{Details}.}
\item{log}{logical, should results be returned in the log scale.}
  \item{scale}{a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm -- see \bold{Details}.}
  \item{weight}{a numeric vector of weights to be applied to the predictions of each subject.}
  \item{init.b}{a numeric matrix of initial values for the random effects. These are used in the optimization procedure that finds the 
mode of the posterior distribution described in Step 2 below.}
  \item{seed}{numeric scalar, the random seed used to produce the results.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  Based on a fitted joint model (represented by \code{object}), and a history of longitudinal responses 
  \eqn{\tilde{y}_i(t) = \{y_i(s), 0 \leq s \leq t\}}{tilde{y_i}(t) = {y_i(s), 0 \leq s \leq t}} and a covariates vector \eqn{x_i} (stored in
  \code{newdata}), this function provides estimates of \eqn{Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)}{Pr(T_i > u | T_i > t, 
  tilde{y}_i(t), x_i)}, i.e., the conditional probability of surviving time \eqn{u} given that subject \eqn{i}, with covariate information
  \eqn{x_i}, has survived up to time \eqn{t} and has provided longitudinal the measurements \eqn{\tilde{y}_i(t)}{tilde{y}_i(t)}.
  
  To estimate \eqn{Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)}{Pr(T_i > u | T_i > t, tilde{y}_i(t), x_i)} and if \code{simulate = TRUE}, a 
  Monte Carlo procedure is followed with the following steps:
  \describe{
    \item{Step 1:}{Take randomly a realization, say \eqn{\theta^*} from the MCMC sample of posterior of the joint model represented by \code{object}.}
    \item{Step 2:}{Simulate random effects values, say \eqn{b_i^*}, from their posterior distribution given survival up to time \eqn{t},
        the vector of longitudinal responses \eqn{\tilde{y}_i(t)} and \eqn{\theta^*}. This is achieved using a Metropolis-Hastings algorithm with
        independent proposals from a properly centered and scaled multivariate \eqn{t} distribution. The \code{scale} argument controls the 
        acceptance rate for this algorithm.}
    \item{Step 3}{Using \eqn{\theta^*} and \eqn{b_i^*}, compute \eqn{Pr(T_i > u | T_i > t, b_i^*, x_i; \theta^*)}{Pr(T_i > 
        u | T_i > t, b_i^*, x_i; \theta^*)}.}
    \item{Step 4:}{Repeat Steps 1-3 \code{M} times.}
  }
  Based on the \code{M} estimates of the conditional probabilities, we compute useful summary statistics, such as their mean, median, and
  percentiles (to produce a confidence interval).

  If \code{simulate = FALSE}, then survival probabilities are estimated using the formula \deqn{Pr(T_i > u | T_i > t, \hat{b}_i, x_i; 
    \hat{\theta}),}{Pr(T_i > u | T_i > t, hat{b}_i, x_i; hat{\theta}),} where \eqn{\hat{\theta}} denotes the posterior means for the parameters, 
    and \eqn{\hat{b}_i} denotes the posterior means for the random effects.
}

\value{
  A list of class \code{survfit.JMbayes} with components:
  \item{summaries}{a list with elements numeric matrices with numeric summaries of the predicted probabilities for each subject.}
  \item{survTimes}{a copy of the \code{survTimes} argument.}
  \item{last.time}{a numeric vector with the time of the last available longitudinal measurement of each subject.}
  \item{obs.times}{a list with elements numeric vectors denoting the timings of the longitudinal measurements for each subject.}
  \item{y}{a list with elements numeric vectors denoting the longitudinal responses for each subject.}
  \item{full.results}{a list with elements numeric matrices with predicted probabilities for each subject in each replication of the Monte Carlo
    scheme described above.}
  \item{success.rate}{a numeric vector with the success rates of the Metropolis-Hastings algorithm described above for each subject.}
  \item{scale}{a copy of the \code{scale} argument.}
}

\references{
Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for longitudinal and time-to-event data. 
\emph{Biometrics} \bold{67}, 819--829.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{plot.survfit.JMbayes}}, \code{\link{predict.JMbayes}}, 
\code{\link{aucJM}}, \code{\link{dynCJM}}, \code{\link{prederrJM}}, \code{\link{jointModelBayes}}}

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

# we will compute survival probabilities for Subject 2 in a dynamic manner, 
# i.e., after each longitudinal measurement is recorded
ND <- pbc2[pbc2$id == 2, ] # the data of Subject 2
survPreds <- vector("list", nrow(ND))
for (i in 1:nrow(ND)) {
    survPreds[[i]] <- survfitJM(jointFit, newdata = ND[1:i, ])
}
survPreds

###########################################################################

# Predictions from multivariate models

pbc2 <- pbc2[!is.na(pbc2$serChol), ]
pbc2.id <- pbc2[!duplicated(pbc2$id), ]
pbc2.id$Time <- pbc2.id$years
pbc2.id$event <- as.numeric(pbc2.id$status != "alive")

# Fit a trivariate joint model
MixedModelFit <- mvglmer(list(log(serBilir) ~ year + (year | id),
                              sqrt(serChol) ~ year + (year | id),
                              hepatomegaly ~ year + (year | id)), data = pbc2,
                         families = list(gaussian, gaussian, binomial), engine = "STAN")

CoxFit <- coxph(Surv(Time, event) ~ drug + age, data = pbc2.id, model = TRUE)

JMFit <- mvJointModelBayes(MixedModelFit, CoxFit, timeVar = "year")

# We want survival probabilities for three subjects
ND <- pbc2[pbc2$id \%in\% c(2, 25, 81), ]

sprobs <- survfitJM(JMFit, ND)
sprobs

# Basic plot
plot(sprobs)

# split in a 2 rows 2 columns and include the survival function in 
# a separate panel; plot only the third & first subjects; change various defaults
plot(sprobs, split = c(2, 2), surv_in_all = FALSE, which_subjects = c(3, 1),
     lty_lines_CI = 3, col_lines = "blue", col_fill_CI = "red", 
     col_points = "pink", pch_points = 12)

###########################################################################

# run Shiny app
runDynPred()
}
}

\keyword{methods}
