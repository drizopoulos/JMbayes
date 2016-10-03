\name{JMbayesObject}

\alias{JMbayesObject}

\title{Fitted JMbayes Object}

\description{
An object returned by the \code{jointModelBayes} function, inheriting from class \code{JMbayes} and representing a fitted
joint model for longitudinal and time-to-event data. Objects of this class have methods for the generic functions 
\code{coef}, \code{confint}, \code{fixed.effects}, \code{logLik}, \code{plot}, \code{print}, 
\code{random.effects}, \code{summary}, and \code{vcov}.
}

\value{
  The following components must be included in a legitimate \code{JMbayes} object.
  \item{mcmc}{a list with the MCMC samples for each parameter (except from the random effects if control argument \code{keepRE} is \code{FALSE}).}
  \item{postMeans}{a list with posterior means.}
  \item{postModes}{a list with posterior modes calculated using kernel desnisty estimation.}
  \item{postVarsRE}{a list with the posterior variance-covariance matrix for the random effects of each subject.}
  \item{StErr}{a list with posterior standard errors.}
  \item{EffectiveSize}{a list with effective sample sizes.}
  \item{StDev}{a list with posterior standard deviations.}
  \item{CIs}{a list with 95\% credible intervals.}
  \item{vcov}{the variance-covariance matrix of the model's parameters based.}
  \item{pD}{the pD value.}
  \item{DIC}{the deviance information criterion value.}
  \item{CPO}{the conditional predictive ordinate value.}
  \item{LPML}{the log pseudo marginal likelihood value.}
  \item{time}{the time used to fit the model.}
  \item{scales}{a list with scaling constants in the Metropolis algorithm.}
  \item{Covs}{a list with the covariance matrices of the proposals in the Metropolis algorithm.}
  \item{acceptRates}{a list of acceptance rates.}
  \item{x}{a list with the design matrices for the longitudinal and event processes.}
  \item{y}{a list with the response vectors for the longitudinal and event processes.}
  \item{Data}{a list of data frames with the data used to fit the models.}
  \item{Terms}{a list of terms objects for the various parts of the joint model.}
  \item{Funs}{a list of functions used for the various parts of the joint model.}
  \item{Forms}{a list of formulas for the two submodels.}
  \item{timeVar}{the value of the \code{timeVar} argument}
  \item{control}{the value of the \code{control} argument.}
  \item{densLongCheck}{a logical indicating whether a scale parameter is required in the longitudinal submodel.}
  \item{param}{the value of the \code{param} argument.}
  \item{priors}{a list with the specification of the prior distributions for the model's parameters. This has the same components as
    the \code{priors} argument of the \code{\link{jointModelBayes}} function.}
  \item{baseHaz}{the value of the \code{baseHaz} argument.}
  \item{df.RE}{the value of the \code{df.RE} argument.}
  \item{call}{the matched call.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{
    \code{\link{jointModelBayes}}
}

\keyword{multivariate}
\keyword{regression}
