\name{logLik.JMbayes}

\alias{logLik.JMbayes}



\title{Log-Likelihood for Joint Models}

\description{
    Computes the log-likelihood for a fitted joint model.
}

\usage{
\method{logLik}{JMbayes}(object, thetas, b, priors = TRUE, marginal.b = TRUE, 
    marginal.thetas = FALSE, full.Laplace = FALSE, useModes = TRUE, \dots)
} 

\arguments{
  \item{object}{an object inheriting from class \code{JMBayes}.}
  \item{thetas}{a list with values for the joint model's parameters. This should have the same structure as
    the \code{coefficients} component of a fitted joint model. If missing \code{object$postMeans} is used.}
  \item{b}{a numeric matrix with random effects value. This should have the same structure as
    the \code{ranef} component of a fitted joint model. If missing \code{ranef(object)} is used.}
  \item{priors}{logical, if \code{TRUE} the priors are also included in the computation.}
  \item{marginal.b}{logical, if \code{TRUE} the marginal log-likelihood over the random effects is returned.
    This marginalization is done using a Laplace approximation.}
  \item{marginal.thetas}{logical, if \code{TRUE} the marginal log-likelihood over the parameters is returned.
    This marginalization is done using a Laplace approximation.}
  \item{full.Laplace}{logical, if \code{FALSE} the posterior means and posterior variances are used in the Laplace
    approximation instead of the posterior modes and posterior hessian matrix of the random effects. Sacrificing a bit of
    accuracy, this will be much faster than calculating the posterior modes. Relevant only when \code{marginal.b = TRUE}.}
  \item{useModes}{logical, if \code{TRUE} the modes are used in the Laplace approximation otherwise the means.}
  \item{\dots}{extra arguments; currently none is used.}
}

\details{
  Let \eqn{y_i} denote the vectors of longitudinal responses, \eqn{T_i} the observed event time, and \eqn{\delta_i} 
  the event indicator for subject \eqn{i} (\eqn{i = 1, \ldots, n}). Let also \eqn{p(y_i | b_i; \theta)} denote the probability
  density function (pdf) for the linear mixed model, \eqn{p(T_i, \delta_i | b_i; \theta)} the pdf for the survival submodel, and
  \eqn{p(b_i; \theta)} the multivariate normal pdf for the random effects, where \eqn{\theta} denotes the full parameter vector. Then,
  if \code{priors = TRUE}, and \code{marginal.b = TRUE}, function \code{logLik()} computes 
  \deqn{\log \int p(y_i | b_i; \theta) p(T_i, \delta_i | b_i; \theta) p(b_i; \theta) db_i + \log p(\theta),}
  where \eqn{p(\theta)} denotes the prior distribution for the parameters. If \code{priors = FALSE} the prior is excluded from the
  computation, i.e., \deqn{\log \int p(y_i | b_i; \theta) p(T_i, \delta_i | b_i; \theta) p(b_i; \theta) db_i,} and when 
  \code{marginal.b = FALSE}, then the conditional on the random effects log-likelihood is computed, i.e.,
   \deqn{\log p(y_i | b_i; \theta) + \log p(T_i, \delta_i | b_i; \theta) + \log p(b_i; \theta) + \log p(\theta),} when 
   \code{priors = TRUE} and \deqn{\log p(y_i | b_i; \theta) + \log p(T_i, \delta_i | b_i; \theta) + \log p(b_i; \theta),}
   when \code{priors = FALSE}.
  
}

\value{
  a numeric scalar of class \code{logLik} with the value of the log-likelihood. It also has
  the attributes \code{df} the number of parameter (excluding the random effects), and \code{nobs}
  the number of subjects. 
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D., Hatfield, L., Carlin, B. and Takkenberg, J. (2014). Combining dynamic predictions from joint models for longitudinal and time-to-event data using Bayesian model averaging. \emph{Journal of the American Statistical Association}. to appear.
}

\seealso{\code{\link{jointModelBayes}}}

\examples{
\dontrun{
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2, 
    random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)

jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

logLik(jointFit)
logLik(jointFit, priors = FALSE)
logLik(jointFit, marginal.b = FALSE)
}
}

\keyword{regression}

