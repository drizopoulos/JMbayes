\name{anova}

\alias{anova.JMbayes}

\title{Anova Method for Fitted Joint Models}

\description{
  Comparison of (non)nested joint models using information criteria.
}

\usage{
\method{anova}{JMbayes}(object, \dots)
}

\arguments{
  \item{object, \dots}{objects inheriting from class \code{JMbayes}.}
}

\value{
  A data frame with rows the different models, and columns the number or parameters in each model, the
  the log pseudo marginal likelihood value, the deviance information criterion value, and the pD value.
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{jointModelBayes}}}

\examples{
\dontrun{
# composite event indicator
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# linear mixed model with natural cubic splines for the time
# effect
lmeFit.pbc1 <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
                   random = ~ ns(year, 2) | id, method = "ML")

# Cox regression model with baseline covariates
coxFit.pbc1 <- coxph(Surv(years, status2) ~ drug * age, data = pbc2.id, x = TRUE)

# the standard joint model fit with only the m_i(t) term in 
# the linear predictor of the survival submodel
jointFit.pbc1 <- jointModelBayes(lmeFit.pbc1, coxFit.pbc1, timeVar = "year")

# we include the time-dependent slopes term
dForm <- list(fixed = ~ 0 + dns(year, 2), random = ~ 0 + dns(year, 2), 
              indFixed = 2:3, indRandom = 2:3)

jointFit.pbc2 <- update(jointFit.pbc1, param = "td-both", extraForm = dForm)

# we include the cumulative effect of the marker
iForm <- list(fixed = ~ 0 + year + ins(year, 2), random = ~ 0 + year + ins(year, 2), 
              indFixed = 1:3, indRandom = 1:3)

jointFit.pbc3 <- update(jointFit.pbc1, param = "td-extra", extraForm = iForm)


# we compare the three models
anova(jointFit.pbc1, jointFit.pbc2, jointFit.pbc3)
}
}

\keyword{methods}
