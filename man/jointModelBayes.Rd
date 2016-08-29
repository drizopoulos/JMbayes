\name{jointModelBayes}
\alias{jointModelBayes}

\title{
Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fits shared parameter joint models for longitudinal and survival outcomes under a Bayesian approach.
}

\usage{
jointModelBayes(lmeObject, survObject, timeVar,  
    param = c("td-value", "td-extra", "td-both", "shared-betasRE", "shared-RE"), 
    extraForm = NULL, baseHaz = c("P-splines", "regression-splines"), 
    transFun = NULL, densLong = NULL, lag = 0, df.RE = NULL, 
    estimateWeightFun = FALSE, weightFun = NULL, init = NULL, 
    priors = NULL, scales = NULL, control = list(), \dots)
}

\arguments{
  \item{lmeObject}{an object of class 'lme' fitted by function \code{lme()} from package \bold{nlme} or by function \code{glmmPQL()} 
        from package \bold{MASS}.}
  \item{survObject}{an object of class 'coxph' fitted by function \code{coxph()} from package \bold{survival}.}
  \item{timeVar}{a character string indicating the time variable in the mixed effects model.}
  \item{param}{a character string specifying the type of association structure between the longitudinal and 
        survival processes. See \bold{Details}.}
  \item{extraForm}{a list with components \code{fixed} a formula representing the fixed-effects part of the user-defined term, 
        \code{indFixed} a numeric vector indicating which fixed effects of \code{lmeObject} are involved in the user-defined term, 
        \code{random} a formula representing the random-effects part of the user-defined term, and \code{indRamdom} a 
        numeric vector indicating which random effects of \code{lmeObject} are involved in the user-defined term. 
        Required only when \code{param = "td-extra"} or \code{param = "td-both"}. See \bold{Examples}.}
  \item{baseHaz}{a character string specifying the type of the baseline hazard function. See \bold{Details}.}
  \item{transFun}{a function or a named list with elements \code{value} and \code{extra} which should be functions. 
        In either case the functions should always have two arguments, namely \code{x} and \code{data} 
        (even when the second one is not needed). The purpose is to transform the \code{value} and/or \code{extra},
        for example including an interaction term, a nonlinear function, etc.}
  \item{densLong}{a function with arguments \code{y}, \code{eta.y}, \code{scale}, \code{log}, and \code{data} that calculates the
        density of the longitudinal outcome. \code{y} denotes the longitudinal responses, \code{eta.y} the linear predictor that
        includes the fixed and random effects, \code{scale} a possible scale parameter (e.g., the measurement error standard deviation),
        \code{log} a logical argument that controls whether the density should be calculated in the log scale, and \code{data} a 
        data frame which may be used to extract variables used in the definition of the density function (e.g., a censoring indicator for
        left censored longitudinal data).}
  \item{lag}{a numeric scalar denoting a lag effect in the time-dependent covariate represented by the mixed model; default is 0.}
  \item{df.RE}{a numeric scalar denoting the number of degrees of freedom for the Student's-\eqn{t} random-effects distribution. If \code{NULL}
        the random effects are assumed to have a multivariate normal distribution.}
  \item{estimateWeightFun}{logical; experimental, not in use yet.}
  \item{weightFun}{a weight function; experimental, not in use yet.}
  \item{init}{a named list of user-specified initial values:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{tau}{the precision parameter from the linear mixed effects model (i.e., \eqn{\tau = 1/\sigma^2} with \eqn{\sigma} 
            denoting the error terms standard deviation).}
        \item{invD}{the inverse variance-covariance matrix of the random effects.}
        \item{b}{a matrix of random effects values.}
        \item{gammas}{the vector of baseline covariates for the survival model.}
        \item{alphas}{the association parameter(s).}
        \item{Dalphas}{the association parameter for the true slopes formulation.}
        \item{Bs.gammas}{the vector of spline coefficients for the spline-approximated baseline risk function.}
    }
    When this list of initial values does not contain some of these components or contains components 
    not of the appropriate length, then the default initial values are used instead.
  }
  \item{priors}{a named list of user-specified prior parameters:
    \describe{
        \item{priorMean.betas}{the prior mean vector of the normal prior for the fixed effects of the linear mixed effects model.}
        \item{priorTau.betas}{the prior precision matrix of the normal prior for the fixed effects of the 
        linear mixed effects model.}
        \item{priorA.tau}{the prior shape parameter of the Gamma prior for the precision parameter of the 
        linear mixed effects model.}
        \item{priorB.tau}{the prior rate parameter of the Gamma prior for the precision parameter of the 
        linear mixed effects model.}
        \item{priorMean.gammas}{the prior mean vector of the normal prior for the regression coefficients of the survival model.}
        \item{priorTau.gammas}{the prior precision matrix of the normal prior for the regression coefficients of 
        the survival model.}
        \item{priorMean.alphas}{the prior mean vector of the normal prior for the association parameter in the survival model.}
        \item{priorTau.alphas}{the prior precision matrix of the normal prior for the association parameter in
        the survival model.}
        \item{priorMean.Dalphas}{the prior mean vector of the normal prior for the slope association parameter 
        in the survival model.}
        \item{priorTau.Dalphas}{the prior precision matrix of the normal prior for the slope association parameter
        in the survival model.}
        \item{priorMean.Bs.gammas}{the prior mean vector of the normal prior for the spline coefficients of the 
        baseline risk function.}
        \item{priorTau.Bs.gammas}{the prior precision matrix of the normal prior for the spline coefficients of 
        the baseline risk function.}
        \item{priorA.tauBs}{the prior shape parameter of the Gamma prior for the precision parameter of the 
        penalty term when \code{baseHaz = "P-splines"}.}
        \item{priorB.tauBs}{the prior rate parameter of the Gamma prior for the precision parameter of the 
        penal term when \code{baseHaz = "P-splines"}.}
        \item{priorR.D}{the prior precision matrix of the Wishart prior for the precision matrix of the random effects.}
        \item{priorK.D}{the degrees of freedom of the Wishart prior for the precision matrix of the random effects.}
    }
  }
  \item{scales}{a named list with names as in \code{init} specifying scaling constants for the proposal distributions 
        in the Metropolis algorithm.}
  \item{control}{a list of control values with components:
    \describe{        
        \item{adapt}{logical default \code{FALSE}; should adaptive Metropolis be used. Currently experimental.}
        \item{n.iter}{integer specifying the total number of iterations; default is 20000.}
        \item{n.burnin}{integer specifying how many of \code{n.iter} to discard as burn-in; default is 3000.}
        \item{n.thin}{integer specifying the thinning of the chains; default is to set \code{n.thin} such that 2000 samples are kept.}
        \item{n.adapt}{integer specifying the number of iterations to use for adaptation; default is 3000.}
        \item{keepRE}{logical; if \code{TRUE} the MCMC samples for the random effect are kept in the output object.}
        \item{priorVar}{integer used as the prior precision in the normal prior for the fixed effects, 
            the regression coefficients of the survival submodel, the association parameters, the extra association parameters, 
            and in the spline coefficients; default is 100.}
        \item{knots}{a numeric vector of knots positions for the spline approximation of the log baseline risk function; 
            default is \code{NULL}, which means that the knots are calculated based on the percentiles of the observed event times.}
        \item{ObsTimes.knots}{logical; if \code{TRUE} (default), the knots are set using the percentiles of the observed event times
            (i.e., including both true events and censored observations). If \code{FALSE}, the knots are set based on the percentiles of
            the true event times alone.}
        \item{lng.in.kn}{a numeric scalar indicating the number of knots to use (based on the percentiles); default is 15 for the penalized
            spline baseline hazard and 5 for the regression spline baseline hazard.}
        \item{ordSpline}{a numeric scalar setting the order of the spline function. This is the number of coefficients in 
            each piecewise polynomial segment, thus a cubic spline has order 4; default is 4.}
        \item{diff}{a numeric scalar setting the order of the differences in the calculation of the penalty term for the penalized 
            baseline hazard; default is 2.}
        \item{seed}{a numeric scalar setting the random seed; default is 1.}
        \item{verbose}{logical; if \code{TRUE} (default), a progress bar is shown in the console.}
    }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{
Function \code{jointModelBayes} fits joint models for longitudinal and survival data under a Bayesian approach. For the longitudinal 
responses a linear mixed effects model represented by the \code{lmeObject} is assumed, unless the user specifies his own probability 
density function using argument \code{densLong}. For the survival times, let \eqn{w_i} denote the vector of baseline covariates in 
\code{survObject}, with associated parameter vector \eqn{\gamma}, \eqn{m_i(t)} the subject-specific linear predictor at time point 
\eqn{t} as defined by the mixed model (i.e., \eqn{m_i(t)} equals the fixed-effects part \code{+} random-effects part of the mixed 
effects model for sample unit \eqn{i}), \eqn{m_i'(t)} denotes an extra user-defined term (based on the specification of argument 
\code{extraForm}) to be included in the linear predictor of the survival submodel, and \eqn{\alpha} and \eqn{\alpha_e} vector of 
association parameters for \eqn{m_i(t)} and \eqn{m_i'(t)}, respectively. Then, \code{jointModelBayes}
assumes a relative risk model of the form \deqn{h_i(t) = h_0(t) \exp\{\gamma^\top w_i + f(m_i(t), m_i'(t), b_i; \alpha, \alpha_e)\},} 
where the baseline risk function is approximated using splines, i.e., \deqn{\log h_0(t) = \sum_k \tilde\gamma_k B(t; \lambda),} 
with \eqn{B(.)} denoting a B-spline basis function, \eqn{\lambda} a vector of knots, and \eqn{\tilde \gamma_k} the corresponding 
splines coefficients (\eqn{\tilde \gamma} correspond to \code{Bs.gammas} above). Argument \code{baseHaz} specifies whether a 
penalized- or regression-spline-approximation is employed. For the former the P-splines approach of Eilers and Marx (1996) is used, 
namely the prior for \eqn{\tilde \gamma} is taken to be proportional to \deqn{p(\tilde \gamma) \propto 
\exp \Bigl(- \frac{\tau_{bs}}{2} \tilde \gamma^\top \Delta^\top \Delta \tilde \gamma \Bigr),} where \eqn{\Delta} denotes the 
differences matrix (the order of the differences is set by the control argument \code{diff}).

Function \eqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d)} specifies the association structure between the two processes. In particular, for
\code{param = "td-value"}, \deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = f_1(m_i(t)) \alpha,} for \code{param = "td-extra"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = f_2(m_i'(t)) \alpha_e,} for \code{param = "td-both"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = f_1(m_i(t)) \alpha + f_2(m_i'(t)) \alpha_e,} for \code{param = "shared-RE"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha^\top b_i,} and for \code{param = "shared-betasRE"},
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha^\top (\beta^* + b_i),} where \eqn{f_1(.)} and \eqn{f_2(.)} denote possible
transformation functions, \eqn{b_i} denotes the vector of random effects for 
the \eqn{i}th subject and \eqn{\beta^*} the fixed effects that correspond to the random effects.
}

\note{
1. The \code{lmeObject} argument should represent a mixed model object without any special structure in the random-effects
covariance matrix (i.e., no use of \code{pdMats}).

2. The \code{lmeObject} object should not contain any within-group correlation structure (i.e., \code{correlation} 
argument of \code{lme()}) or within-group heteroscedasticity structure (i.e., \code{weights} argument of \code{lme()}).

3. It is assumed that the linear mixed effects model \code{lmeObject} and the survival model \code{survObject} have been 
fitted to the same subjects. Moreover, it is assumed that the ordering of the subjects is the same for both 
\code{lmeObject} and \code{survObject}, i.e., that the first line in the data frame containing the event times 
corresponds to the first set of lines identified by the grouping variable in the data frame containing the repeated 
measurements, and so on. Furthermore, the scale of the time variable (e.g., days, months, years) should be the same in both
the \code{lmeObject} and \code{survObject} objects.

4. In the \code{print} and \code{summary} generic functions for class \code{jointModel}, the estimated coefficients (and
standard errors for the \code{summary} generic) for the event process are augmented with the element "Assoct" that 
corresponds to the association parameter \eqn{\alpha} and the element "AssoctE" that corresponds to the parameter 
\eqn{\alpha_e} when \code{parameterization} is \code{"td-extra"} or \code{"td-both"} (see \bold{Details}).
}

\value{
See \code{\link{JMbayesObject}} for the components of the fit.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements and event time data.
\emph{Biostatistics} \bold{1}, 465--480.

Hsieh, F., Tseng, Y.-K. and Wang, J.-L. (2006) Joint modeling of survival and longitudinal data: Likelihood 
approach revisited. \emph{Biometrics} \bold{62}, 1037--1043.

Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: With 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011) Dynamic predictions and prospective accuracy in joint models for longitudinal 
and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\seealso{
\code{\link{coef.JMbayes}},
\code{\link{ranef.JMbayes}}, 
\code{\link{logLik.JMbayes}}, 
\code{\link{survfitJM}},
\code{\link{aucJM}},
\code{\link{dynCJM}},
\code{\link{prederrJM}}, 
\code{\link{predict.JMbayes}}
}

\examples{

\dontrun{
# A joint model for the AIDS dataset:
# First we fit the linear mixed model for the longitudinal measurements of
# sqrt CD4 cell counts
lmeFit.aids <- lme(CD4 ~ obstime * drug, random = ~ obstime | patient, data = aids)
# next we fit the Cox model for the time to death (including the 'x = TRUE' argument)
survFit.aids <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# the corresponding joint model is fitted by (the default is to assume 
# the current value parameterization)
jointFit.aids <- jointModelBayes(lmeFit.aids, survFit.aids, timeVar = "obstime")
summary(jointFit.aids)

# A joint model for the PBC dataset:
# We first fit the linear mixed and Cox models. In the first we include 
# splines to model flexibly the subject-specific longitudinal trajectories
lmeFit.pbc <- lme(log(serBilir) ~ ns(year, 2),
    random = list(id = pdDiag(form = ~ ns(year, 2))), data = pbc2)
survFit.pbc <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)

# the corresponding joint model is fitted by:
jointFit.pbc <- jointModelBayes(lmeFit.pbc, survFit.pbc, timeVar = "year", 
    baseHaz = "regression-splines")
summary(jointFit.pbc)

# we update the joint model fitted for the PBC dataset by including
# the time-dependent slopes term. To achieve this we need to define 
# the 'extraForm' argument, in which we use function dns() to numerically
# compute the derivative of the natural cubic spline. In addition, we increase
# the number of MCMC iterations to 35000
dform = list(fixed = ~ 0 + dns(year, 2), random = ~ 0 + dns(year, 2),
    indFixed = 2:3, indRandom = 2:3)
jointFit.pbc2 <- update(jointFit.pbc, param = "td-both", extraForm = dform,
    n.iter = 35000)
summary(jointFit.pbc2)

# we fit the same model with the shared random effects formulation
jointFit.pbc3 <- update(jointFit.pbc, param = "shared-betasRE")
summary(jointFit.pbc3)

# a joint model for left censored longitudinal data
# we create artificial left censoring in serum bilirubin
pbc2$CensInd <- as.numeric(pbc2$serBilir <= 0.8)
pbc2$serBilir2 <- pbc2$serBilir
pbc2$serBilir2[pbc2$CensInd == 1] <- 0.8

censdLong <- function (y, eta.y, scale, log = FALSE, data) {
    log.f <- dnorm(x = y, mean = eta.y, sd = scale, log = TRUE)
    log.F <- pnorm(q = y, mean = eta.y, sd = scale, log.p = TRUE)
    ind <- data$CensInd
    if (log) {
        (1 - ind) * log.f + ind * log.F
    } else {
        exp((1 - ind) * log.f + ind * log.F)
    }
}
lmeFit.pbc2 <- lme(log(serBilir2) ~ ns(year, 2), data = pbc2,
                   random = ~ ns(year, 2) | id, method = "ML")
jointFit.pbc4 <- jointModelBayes(lmeFit.pbc2, survFit.pbc, timeVar = "year",
                                  densLong = censdLong)

summary(jointFit.pbc4)

# a joint model for a binary outcome
pbc2$serBilirD <- 1 * (pbc2$serBilir > 1.8)
lmeFit.pbc3 <- glmmPQL(serBilirD ~ year, random = ~ year | id, 
	family = binomial, data = pbc2)

dLongBin <- function (y, eta.y, scale, log = FALSE, data) {
    dbinom(x = y, size = 1, prob = plogis(eta.y), log = log)
}

jointFit.pbc5 <- jointModelBayes(lmeFit.pbc3, survFit.pbc, timeVar = "year", 
	densLong = dLongBin)

summary(jointFit.pbc5)


# create start-stop counting process variables
pbc <- pbc2[c("id", "serBilir", "drug", "year", "years",
              "status2", "spiders")]
pbc$start <- pbc$year
splitID <- split(pbc[c("start", "years")], pbc$id)
pbc$stop <- unlist(lapply(splitID,
                          function (d) c(d$start[-1], d$years[1]) ))
pbc$event <- with(pbc, ave(status2, id,
                           FUN = function (x) c(rep(0, length(x)-1), x[1])))
pbc <- pbc[!is.na(pbc$spiders), ]

# left-truncation
pbc <- pbc[pbc$start != 0, ] 

lmeFit.pbc <- lme(log(serBilir) ~ drug * ns(year, 2),
                  random = ~ ns(year, 2) | id, data = pbc)

tdCox.pbc <- coxph(Surv(start, stop, event) ~ drug * spiders + cluster(id),
                   data = pbc, x = TRUE, model = TRUE)

jointFit.pbc6 <- jointModelBayes(lmeFit.pbc, tdCox.pbc, timeVar = "year")

summary(jointFit.pbc6)
}
}

\keyword{multivariate}
\keyword{regression}
