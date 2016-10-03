\name{mvJointModelBayes}
\alias{mvJointModelBayes}

\title{
Multivariate Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fits multivariate shared parameter joint models for longitudinal and survival outcomes under a Bayesian approach.
}

\usage{
mvJointModelBayes(mvglmerObject, coxphObject, timeVar,
    Formulas = list(NULL), Interactions = list(NULL),
    priors = NULL, control = NULL, \dots)
}

\arguments{
  \item{mvglmerObject}{an object of class 'mvglmer' fitted by function 
    \code{mvglmer()}.}
  
  \item{coxphObject}{an object of class 'coxph' fitted by function \code{coxph()} 
    from package \bold{survival}.}
  
  \item{timeVar}{a character string indicating the time variable in the multivariate 
    mixed effects model.}
  
  \item{Formulas}{a list of lists. Each inner list should have components 
    \code{fixed} a formula representing the fixed-effects part of the user-defined term,
    \code{indFixed} a numeric vector indicating which fixed effects of 
    \code{mvglmerObject} are involved in the user-defined term, 
    \code{random} a formula representing the random-effects part of the user-defined term, 
    and \code{indRamdom} a numeric vector indicating which random effects of 
    \code{mvglmerObject} are involved in the user-defined term. See \bold{Examples}.}
  
  \item{Interactions}{a list specifying interaction terms for the components of the 
    longitudinal outcomes that are included in the survival submodel. See \bold{Examples}.}
  
  \item{priors}{a named list of user-specified prior parameters:
    \describe{
        \item{mean_Bs_gammas}{the prior mean vector of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{Tau_Bs_gammas}{the prior precision matrix of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{mean_gammas}{the prior mean vector of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{Tau_gammas}{the prior precision matrix of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{mean_alphas}{the prior mean vector of the normal prior for the association
            parameters.}
        \item{Tau_alphas}{the prior mean vector of the normal prior for the association
            parameters.}
        \item{A_tau_Bs_gammas}{the prior shape parameter of the Gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{B_tau_Bs_gammas}{the prior rate parameter of the Gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{shrink_gammas}{logical; should the regression coefficients for the baseline
            covariates be shrinked.}
        \item{A_tau_gammas}{the prior shape parameter of the Gamma prior for the
            precision parameter of the global penalty term baseline regression coefficients.
            Only relevant when \code{shrink_gammas = TRUE}.}
        \item{B_tau_gammas}{the prior rate parameter of the Gamma prior for the
            precision parameter of the penalty term for the baseline regression coefficients. 
            Only relevant when \code{shrink_gammas = TRUE}.}
        \item{A_phi_gammas}{the prior shape parameter of the Gamma prior for the
            precision parameters of each baseline regression coefficients.
            Only relevant when \code{shrink_gammas = TRUE}.}
        \item{B_phi_gammas}{the prior rate parameter of the Gamma prior for the
            precision parameters of each baseline regression coefficients.
            Only relevant when \code{shrink_gammas = TRUE}.}
        \item{shrink_alphas}{logical; should the association parameters be shrinked.}
        \item{A_tau_alphas}{the prior shape parameter of the Gamma prior for the
            precision parameter of the global penalty term for the association parameters.
            Only relevant when \code{shrink_alphas = TRUE}.}
        \item{B_tau_alphas}{the prior rate parameter of the Gamma prior for the
            precision parameter of the penalty term for the association parameters. 
            Only relevant when \code{shrink_alphas = TRUE}.}
        \item{A_phi_alphas}{the prior shape parameter of the Gamma prior for the
            precision parameters of each association parameter.
            Only relevant when \code{shrink_alphas = TRUE}.}
        \item{B_phi_alphas}{the prior rate parameter of the Gamma prior for the
            precision parameters of each association parameter.
            Only relevant when \code{shrink_alphas = TRUE}.}
    }
  }
  
  \item{control}{a list of control values with components:
    \describe{        
        \item{n_iter}{integer specifying the total number of iterations after burn in; 
            default is 300.}
        \item{n_burnin}{integer specifying how many of iterations to discard as burn-in; 
            default is 1000.}
        \item{n_thin}{integer specifying the thinning of the chains; 
            default is 300.}
        \item{n_block}{integer specifying the number of block iterations in which the 
            acceptance rates are checked; default is 100.}
        \item{target_acc}{a numeric scalar denoting the target acceptance rate; 
            default is 0.234.}
        \item{knots}{a numeric vector of knots positions for the spline approximation of 
            the log baseline risk function; default is \code{NULL}, which means that the 
            knots are calculated based on the percentiles of the observed event times.}
        \item{ObsTimes.knots}{logical; if \code{TRUE} (default), the knots are set using 
            the percentiles of the observed event times (i.e., including both true events 
            and censored observations). If \code{FALSE}, the knots are set based on the 
            percentiles of the true event times alone.}
        \item{lng.in.kn}{a numeric scalar indicating the number of knots to use (based on 
            the percentiles); default is 15.}
        \item{ordSpline}{an integer setting the order of the spline function. This is the 
            number of coefficients in each piecewise polynomial segment, thus a cubic 
            spline has order 4; default is 4.}
        \item{diff}{an integer setting the order of the differences in the calculation of 
            the penalty term for the penalized baseline hazard; default is 2.}
        \item{seed}{an integer setting the random seed; default is 1.}
        \item{n_cores}{an integer specifying the number of cores to use. Default is the
            the available cores minus one.}
        \item{GQsurv}{a character string specifying the type of Gaussian quadrature to 
            be used. Options are "GaussKronrod" (default) and "GaussLegendre".}
        \item{GQsurv.k}{an integer specifying the number of quadrature points; 
            default is 15.}
    }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{

** to fill in ** 
}

\value{
A list of class \code{mvJMbayes} with components:

\item{mcmc}{a list with the MCMC samples for each parameter.}
\item{components}{a copy of the \code{components} element of \code{mvglmerObject}.}
\item{Data}{a list of data used to fit the model.}
\item{control}{a copy of the \code{control} values used in the fit.}
\item{mcmc.info}{a list with information over the MCMC (i.e., time it took, iterations, 
    etc.).}
\item{priors}{a copy of the priors used.}
\item{postwMeans}{a list with posterior weighted means.}
\item{postMeans}{a list with posterior means.}
\item{postModes}{a list with posterior modes calculated using kernel desnisty estimation.}
\item{EffectiveSize}{a list with effective sample sizes.}
\item{StErr}{a list with posterior standard errors.}
\item{StDev}{a list with posterior standard deviations.}
\item{CIs}{a list with 95\% credible intervals.}
\item{Pvalues}{a list of tail probabilities for containg the zero value.}
\item{call}{the matched call.}
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
\code{\link{mvglmer}},
\code{\link{jointModelBayes}}
}

\examples{

\dontrun{
pbc2.id$Time <- pbc2.id$years
pbc2.id$event <- as.numeric(pbc2.id$status != "alive")

##########################################################################################

##############################################
# Univariate joint model for serum bilirubin #
##############################################

# [1] Fit the mixed model using mvglmer(). The main arguments of this function are:
# 'formulas' a list of lme4-like formulas (a formular per outcome),
# 'data' a data.frame that contains all the variables specified in 'formulas' (NAs allowed),
# 'families' a list of family objects specifying the type of each outcome (currently only
# gaussian, binomial and poisson are allowed).
MixedModelFit1 <- mvglmer(list(log(serBilir) ~ year + (year | id)), data = pbc2,
                          families = list(gaussian))

# [2] Fit a Cox model, specifying the baseline covariates to be included in the joint
# model; you need to set argument 'model' to TRUE.
CoxFit1 <- coxph(Surv(Time, event) ~ drug + age, data = pbc2.id, model = TRUE)

# [3] The basic joint model is fitted using a call to mvJointModelBayes(), which is very
# similar to JointModelBayes(), i.e.,
JMFit1 <- mvJointModelBayes(MixedModelFit1, CoxFit1, timeVar = "year")
summary(JMFit1)
plot(JMFit1)

##########################################################################################

#########################################################
# Bivariate joint model for serum bilirubin and spiders #
#########################################################

MixedModelFit2 <- mvglmer(list(log(serBilir) ~ year + (year | id),
                               spiders ~ year + (1 | id)), data = pbc2,
                          families = list(gaussian, binomial))

CoxFit2 <- coxph(Surv(Time, event) ~ drug + age, data = pbc2.id, model = TRUE)

JMFit2 <- mvJointModelBayes(MixedModelFit2, CoxFit2, timeVar = "year")
summary(JMFit2)
plot(JMFit2)

##########################################################################################

#######################
# slopes & area terms #
#######################

# We extend model 'JMFit2' by including the value and slope term for bilirubin, and
# the area term for spiders (in the log-odds scale). To include these terms into the model
# we specify the 'Formulas' argument. This is specified in a similar manner as the
# 'derivForms' argument of jointModelBayes(). In particular, it should be a list of lists.
# Each component of the outer list should have as name the name of the corresponding
# outcome variable. Then in the inner list we need specify four components, namely,
# 'fixed' & 'random' R formulas specifying the fixed and random effects part of the term
# to be included, and 'indFixed' & 'indRandom' integer indicices specifying which of the
# original fixed and random effects are involved in the claculation of the new term. In
# the inner list you can also optionally specify a name for the term you want to include.
# Notes: (1) For terms not specified in the 'Formulas' list, the default value functional
# form is used. (2) If for a particular outcome you want to include both the value
# functional form and an extra term, then you need to specify that in the 'Formulas'
# using two entries. To include the value functional form you only need to set the
# corresponding to 'value', and for the second term to specify the inner list. See
# example below on how to include the value and slope for serum bilirubin (for example,
# if the list below the entry '"log(serBilir)" = "value"' was not give, then only the
# slope term would have been included in the survival submodel).

Forms <- list("log(serBilir)" = "value",
              "log(serBilir)" = list(fixed = ~ 1, random = ~ 1,
                                     indFixed = 2, indRandom = 2, name = "slope"),
              "spiders" = list(fixed = ~ 0 + year + I(year^2/2), random = ~ 0 + year,
                               indFixed = 1:2, indRandom = 1, name = "area"))

JMFit3 <- update(JMFit2, Formulas = Forms)
summary(JMFit3)
plot(JMFit3)

##########################################################################################

#####################
# Interaction terms #
#####################

# We further extend the previous model by including the interactions terms between the
# terms specified in 'Formulas' and 'drug'. The names specified in the list that defined
# the interaction factors should match with the names of the output from 'JMFit3'; the
# only exception is with regard to the 'value' functional form. See specification below
# (to include the interaction of the value term of 'log(serBilir)' with 'drug', in the
# list we can either specify as name of the corresponding formula 'log(serBilir)_value'
# or just 'log(serBilir)'):

Ints <- list("log(serBilir)" = ~ drug, "log(serBilir)_slope" = ~ drug,
             "spiders_area" = ~ drug)

# because of the many association parameters we have, we place a shrinkage prior on the
# alpha coefficients. In particular, if we have K association parameters, we assume that
# alpha_k ~ N(0, tau * phi_k), k = 1, ..., K. The precision parameters tau and phi_k are
# given Gamma priors. Precision tau is global shrinkage parameter, and phi_k a specific
# per alpha coefficient.
JMFit4 <- update(JMFit3, Interactions = Ints, priors = list(shrink_alphas = TRUE))
summary(JMFit4)
plot(JMFit4)

}
}

\keyword{multivariate}
\keyword{regression}
