\name{JMbayes}
\alias{JMbayes-package}
\alias{JMbayes}
\docType{package}

\title{
Joint Modeling of Longitudinal and Time-to-Event Data in R under a Bayesian Approach
}

\description{
This package fits shared parameter models for the joint modeling of normal longitudinal responses and event times 
under a Bayesian approach. Various options for the survival model and the association structure are
provided.
}

\details{
\tabular{ll}{
Package: \tab JMbayes\cr
Type: \tab Package\cr
Version: \tab 0.7-6\cr
Date: \tab 2015-09-01\cr
License: \tab GPL (>=2)\cr
}

The package has a single model-fitting function called \code{\link{jointModelBayes}}, which accepts as main arguments a linear
mixed effects object fit returned by function \code{lme()} of package \pkg{nlme}, and a Cox model object fit returned 
by function \code{coxph()} of package \pkg{survival}. The \code{survMod} argument of specifies the type of survival submodel 
to be fitted; available options are a relative risk model with a Weibull baseline hazard (default) and a relative risk model
with a B-spline approximation of the log baseline risk function. In addition, the \code{param} specifies the association structure
between the longitudinal and survival processes; available options are: \code{"td-value"} which is the classic formulation used in 
Wulfsohn and Tsiatis (1997); \code{"td-extra"} which is a user-defined, possibly time-dependent, term based on the specification of 
the \code{extraForm} argument of \code{\link{jointModelBayes}}. This could be used to include terms, such as the time-dependent 
slope (i.e., the derivative of the subject-specific linear predictor of the linear mixed model) and the time-dependent cumulative 
effect (i.e., the integral of the subject-specific linear predictor of the linear mixed model); \code{"td-both"} which is the 
combination of the previous two parameterizations, i.e., the current value and the user-specified terms are included in the linear 
predictor of the relative risk model; and \code{"shared-RE"} where only the random effects of the linear mixed model are included 
in the linear predictor of the survival submodel.

The package also offers several utility functions that can extract useful information from fitted joint models. The most
important of those are included in the \bold{See also} Section below.
}

\author{
Dimitris Rizopoulos

Maintainer: Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>
}

\references{
Guo, X. and Carlin, B. (2004) Separate and joint modeling of longitudinal and event time data 
using standard computer packages. \emph{The American Statistician} \bold{54}, 16--24.

Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements 
and event time data. \emph{Biostatistics} \bold{1}, 465--480.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011) Dynamic predictions and prospective accuracy in joint models for longitudinal 
and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D. and Ghosh, P. (2011) A Bayesian semiparametric multivariate joint model for multiple 
longitudinal outcomes and a time-to-event. \emph{Statistics in Medicine} \bold{30}, 1366--1380.

Rizopoulos, D., Verbeke, G. and Molenberghs, G. (2010) Multiple-imputation-based residuals and diagnostic 
plots for joint models of longitudinal and survival outcomes. \emph{Biometrics} \bold{66}, 20--29.

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\keyword{ package }
\keyword{ multivariate }

\seealso{
    \code{\link{jointModelBayes}},
    \code{\link{survfitJM}},
    \code{\link{aucJM}},
    \code{\link{dynCJM}},
    \code{\link{prederrJM}},
    \code{\link{predict.JMbayes}},
    \code{\link{logLik.JMbayes}}
}

