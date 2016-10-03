\name{mvglmer}
\alias{mvglmer}

\title{
Multivariate Mixed Models
}

\description{
Fits multivariate mixed models under a Bayesian approach using JAGS.
}

\usage{
mvglmer(formulas, data, families, overdispersion = FALSE,
    priors = NULL, control = NULL, \dots)
}

\arguments{
  \item{formulas}{a list of R formulas representing the mixed models; 
        these should be \bold{lme4}-type formulas.}
  
  \item{data}{a data.frame that contains all the variable to be used when fitting the
        multivariate mixed model.}
  
  \item{families}{a list of families objects correspond to each outcome.}
  
  \item{overdispersion}{logical; for Poisson outcomes, should an overdispersion parameter
        be included.}

  \item{priors}{a named list of user-specified prior parameters:
    \describe{
        \item{taus_betas}{the prior precision parameter for the fixed effects; default
            is 0.001.}
        
        \item{priorK.D}{degrees of freedom for the wishart prior for the inverse
            covariance matrix of the random effects; default is number of random effects 
            plus one.}
            
        \item{priorR.D}{precision matrix of the wishart prior for the inverse
            covariance matrix of the random effects; default to a diagonal matrix with
            diagonal ellements given a Gamma prior with parameters \code{A_R.D} and 
            \code{A_R.D}.}
        \item{A_R.D}{the prior shape parameter of the Gamma prior for the diagonal 
            elements of the precision matrix of the wishart prior for the inverse 
            covariance matrix of the random effects; default is 0.5.}
        \item{B_R.D}{the prior shape parameter of the Gamma prior for the diagonal 
            elements of the precision matrix of the wishart prior for the inverse 
            covariance matrix of the random effects; default is 0.001.}
        \item{tau_half_cauchy}{prior precision parameter of a half-Cauchy distribution for 
            the precision parameter of a random intercept, when only a single outcome is
            specified with a single random effect; default is 0.1.}
            
        \item{A_tau}{the prior shape parameter for the precision of the error terms of
            Gaussian outcomes.}
        \item{B_tau}{the prior rate parameter for the precision of the error terms of
            Gaussian outcomes.}
    }
  }
  
  \item{control}{a list of control values with components:
    \describe{        
        \item{n.iter}{integer specifying the total number of iterations after burn in; 
            default is 28000.}
        \item{n.burnin}{integer specifying how many of iterations to discard as burn-in; 
            default is 3000.}
        \item{n.thin}{integer specifying the thinning of the chains; 
            default is 50.}
        \item{n.adapt}{integer specifying the number of adapt iterations in which the 
            acceptance rates are checked; default is 3000.}
        \item{n.chains}{integer specifying the number of chains to use; default is 2.}
        \item{n.processors}{integer specifying the number of processors to use; 
            default is the number of available processors minus one.}
        \item{working.directory}{a character string giving the path on where to save
            the JAGS model; default is the working directory.}
        \item{clear.model}{logical; should the JAGS models be deleted after the model has
            run; default is \code{TRUE}.}
        
        \item{seed}{an integer setting the random seed; default is 1.}
            
     }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{
This function creates a JAGS program representing a multivariate mixed effects that is run
with JAGS using the \bold{jagsUI} package. Currently only Gaussian, Bernoulli and Poisson
longitudinal outcomes can be handled. 
}

\value{
A list of class \code{mvglmer} with components:

\item{mcmc}{a list with the MCMC samples for each parameter.}
\item{components}{a list with design matrices and responses vectors extracted by applying
    the \code{formulas} in \code{data}.}
\item{data}{a copy of \code{data}.}
\item{control}{a copy of the \code{control} values used in the fit.}
\item{mcmc.info}{a list with information over the MCMC (i.e., time it took, iterations, 
    etc.).}
\item{DIC}{the DIC value for the fitted model.}
\item{pD}{the pD valu for the fitted model.}
\item{Rhat}{a list with the Rhat convergence diagnostics for each parameter.}
\item{priors}{a copy of the priors used.}
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
***
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\seealso{
\code{\link{mvJointModelBayes}},
\code{\link{jointModelBayes}}
}

\examples{

\dontrun{
MixedModelFit <- mvglmer(list(log(serBilir) ~ year + (year | id),
                               spiders ~ year + (1 | id)), data = pbc2,
                          families = list(gaussian, binomial))

summary(MixedModelFit)
plot(MixedModelFit)
}
}

\keyword{multivariate}
\keyword{regression}
