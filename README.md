JMbayes: Joint Models for Longitudinal and Survival Data under the Bayesian Approach
================
[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/JMbayes.svg?branch=master)](https://travis-ci.org/drizopoulos/JMbayes) [![CRAN status](http://www.r-pkg.org/badges/version/JMbayes)](https://cran.r-project.org/package=JMbayes) [![](https://cranlogs.r-pkg.org/badges/grand-total/JMbayes)](https://CRAN.R-project.org/package=JMbayes) [![Download counter](http://cranlogs.r-pkg.org/badges/JMbayes)](https://cran.r-project.org/package=JMbayes)
[![Research software impact](http://depsy.org/api/package/cran/JMbayes/badge.svg)](http://depsy.org/package/r/JMbayes)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1038971.svg)](https://doi.org/10.5281/zenodo.1038971)

Description
------------

This repository contains the source files for the R package <strong>JMbayes</strong>. 
This package fits joint models for longitudinal and time-to-event data under a Bayesian 
approach using MCMC. These models are applicable in mainly two settings. First, when focus
is on the survival outcome and we wish to account for the effect of an endogenous 
(aka internal) time-dependent covariates measured with error. Second, when focus is on the
longitudinal outcome and we wish to correct for nonrandom dropout.

The package contains two main joint-model-fitting functions, `jointModelBayes()` and 
`mvJointModelBayes()` with similar syntax but different capabilities.

Basic Features `jointModelBayes()`
------------

- It can fit joint models for a single longitudinal outcome and a time-to-event outcome. 

- The user can specify her own density function for the longitudinal responses using 
argument `densLong` (default is the normal pdf). Among others, this allows to fit joint 
models with categorical and left-censored longitudinal responses and robust joint models 
with Student's-t error terms. In addition, using the `df.RE` argument, the user can also 
change the distribution of the random effects from multivariate normal to a multivariate 
Student's-t with prespecified degrees of freedom.

- For the survival outcome a relative risk models is assumed with a B-spline approximation
for the baseline hazard (penalized (default) or regression splines can be used). 
Left-truncation and exogenous time-varying covariates can also be accommodated.

- The user has now the option to define custom transformation functions for the terms of 
the longitudinal submodel that enter into the linear predictor of the survival submodel 
(arguments `extraForm`, `param`). For example, the current value of the 
longitudinal outcomes, the velocity of the longitudinal outcome (slope), the area under
the longitudinal profile. From the aforementioned options, in each model up to two terms 
can be included. In addition, using argument `transFun` interactions terms, nonlinear terms 
(polynomials, splines) can be considered.

Basic Features `mvJointModelBayes()`
------------

- It can fit joint models for multiple longitudinal outcomes and a time-to-event outcome. 

- The longitudinal part of the joint model is a multivariate generalized linear mixed 
effects models, currently allowing for normal, binary and Poisson outcomes. This model is
first fitted using function `mvglmer()`.

- For the survival outcome a relative risk models is assumed with a B-spline approximation
for the baseline hazard (penalized (default) or regression splines can be used). 
Left-truncation, interval censored data and exogenous time-varying covariates can also be 
accommodated.

- The user has now the option to define custom transformation functions for the terms of 
the longitudinal submodel that enter into the linear predictor of the survival submodel 
(argument `Formulas`). For example, the current value of the longitudinal outcomes, the 
velocity of the longitudinal outcome (slope), the area under the longitudinal profile. 
From the aforementioned options, in each model limitless terms can be included. In 
addition, using argument `Interactions` allows to include interactions terms of the 
longitudinal components with other observed factors. A special case for this argument is
to use function `tve()` that allows for time-varying regression coefficients in the 
relative risk model. Furthermore, argument `transFuns` allows to transform the longitudinal 
components using some pre-defined transformation function (i.e., `exp()`, `expit()`, `log`,
`sqrt()`).

- The aforementioned features are illustrated in the [Multivariate Joint Models vignette](http://www.drizopoulos.com/vignettes/Multivariate%20Joint%20Models.html).

Dynamic predictions
------------

* Function `survfitJM()` computes dynamic survival probabilities.

* Function `predict()` computes dynamic predictions for the longitudinal outcome.

* Function `aucJM()` calculates time-dependent AUCs for joint models, and function 
`rocJM()` calculates the corresponding time-dependent sensitivities and specifies.

* Function `prederrJM()` calculates prediction errors for joint models.

* Function `runDynPred()` invokes a [shiny](https://shiny.rstudio.com/) application that 
can be used to streamline the calculation of dynamic predictions for models fitted by
**JMbayes**.

Vignettes
------------
Vignettes are available in the `doc` directory:

* [Multivariate_Joint_Models.html](http://www.drizopoulos.com/vignettes/multivariate%20joint%20models) illustrates the
basic capabilities of `mvJointModelBayes()`.

* [Dynamic_Predictions.html](http://www.drizopoulos.com/vignettes/dynamic_predictions) illustrates how dynamic 
predictions from multivariate joint models can be computed and evaluated.


