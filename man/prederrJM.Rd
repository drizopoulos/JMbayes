\name{prederrJM}

\alias{prederrJM}
\alias{prederrJM.JMbayes}

\title{Prediction Errors for Joint Models}

\description{ 
  Using the available longitudinal information up to a starting time point, this function computes an estimate 
  of the prediction error of survival at a horizon time point based on joint models.
}

\usage{
prederrJM(object, newdata, Tstart, Thoriz, \dots)

\method{prederrJM}{JMbayes}(object, newdata, Tstart, Thoriz, 
    lossFun = c("square", "absolute"), interval = FALSE, idVar = "id", 
    simulate = FALSE, M = 100, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{newdata}{a data frame that contains the longitudinal and covariate information for the subjects for which prediction 
    of survival probabilities is required. The names of the variables in this data frame must be the same as in the data frames that
    were used to fit the linear mixed effects model (using \code{lme()}) and the survival model (using \code{coxph()})
    that were supplied as the two first argument of \code{\link{jointModelBayes}}. In addition, this data frame should contain a variable 
    that identifies the different subjects (see also argument \code{idVar}).}
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} mast be later than \code{Tstart}.}
  \item{lossFun}{either the options \code{"absolute"} or \code{"square"} (default),  or a user-specified loss function. 
    As the names suggest, 
    when \code{lossFun = "absolute"} the loss function is \eqn{L(x) = |x|}, whereas when \code{lossFun = "square"} the loss function is 
    \eqn{L(x) = x^2}. If a user-specified function is supplied, this should have a single argument and be vectorized.}
  \item{interval}{logical; if \code{TRUE} the weighted prediction error in the interval \code{[Tstart, Thoriz]} is calculated, while
    if \code{FALSE} the prediction error at time \code{Thoriz} is calculated using the longitudinal information up to time \code{Tstart}.}  
  \item{idVar}{the name of the variable in \code{newdata} that identifies the different subjects.}
  \item{simulate}{logical; if \code{TRUE}, a Monte Carlo approach is used to estimate survival probabilities. If \code{FALSE},
    a first order estimator is used instead. See \code{\link{survfitJM}} for mote details.}
  \item{M}{a numeric scalar denoting the number of Monte Carlo samples; see \code{\link{survfitJM}} for mote details.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  Based on a fitted joint model (represented by \code{object}) and using the data supplied in argument \code{newdata}, this function
  computes the following estimate of the prediction: \deqn{PE(u | t) = \{R(t)\}^{-1} \sum_{i: T_i \geq s} I(T_i \geq u) 
  L\{1 - Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)\}}
  \deqn{+ \delta_i I(T_i < u) L\{0 - Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)\}} 
  \deqn{+ (1 - \delta_i) I(T_i < u) [S_i(u \mid T_i, \tilde{y}_i(t)) L\{1 - Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)\}}
  \deqn{+ \{1 - S_i(u \mid T_i, \tilde{y}_i(t))\} L\{0 - Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)\}],}
  where \eqn{R(t)} denotes the number of subjects at risk at time \eqn{t = } \code{Tstart}, \eqn{\tilde{y}_i(t) = \{y_i(s), 0 \leq s \leq t\}} denotes the available 
  longitudinal measurements up to time \eqn{t}, \eqn{T_i} denotes the observed event time for subject \eqn{i}, \eqn{\delta_i} is the event indicator, 
  \eqn{s} is the starting time point \code{Tstart} up to which the longitudinal information is used, and \eqn{u > s} is the horizon time point \code{Thoriz}.
  Function \eqn{L(.)} is the loss function that can be the absolute value (i.e., \eqn{L(x) = |x|}), the squared value (i.e., \eqn{L(x) = x^2}), 
  or a user-specified function. The probabilities \eqn{Pr(T_i > u | T_i > t, \tilde{y}_i(t), x_i)} are calculated by \code{\link{survfitJM}}.
  
  When \code{interval} is set to \code{TRUE}, then function \code{prederrJM} computes the integrated prediction error in the interval 
  \eqn{(u,t) =} \code{(Tstart, Thoriz)} defined as \deqn{IPE(u | t) = \sum_{i: t \leq T_i \leq u} w_i(T_i) PE(T_i | t),} where 
  \deqn{w_i(T_i) = \frac{\delta_i G(T_i) / G(t)}{\sum_{i: t \leq T_i \leq u} \delta_i G(T_i) / G(t)},} with \eqn{G(.)} denoting
  the Kaplan-Meier estimator of the censoring time distribution.
}

\value{
  A list of class \code{prederrJM} with components:
  \item{prederr}{a numeric scalar denoting the estimated prediction error.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}  
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{interval}{a copy of the \code{interval} argument.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}
  \item{lossFun}{a copy of the \code{lossFun} argument.}
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2002). Identification and efficacy of 
longitudinal markers for survival. \emph{Biostatistics} \bold{3}, 33--50.

Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for 
longitudinal and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{survfitJM}}, \code{\link{aucJM}}, \code{\link{dynCJM}}, \code{\link{jointModelBayes}}}

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

# prediction error at year 10 using longitudinal data up to year 5 
prederrJM(jointFit, pbc2, Tstart = 5, Thoriz = 10)
prederrJM(jointFit, pbc2, Tstart = 5, Thoriz = 6.5, interval = TRUE)
}
}

\keyword{methods}
