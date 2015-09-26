# JMbayes
This repository contains the source files for the R package <strong>JMbayes</strong>. 
This package fits joint models for longitudinal and time-to-event data under a Bayesian 
approach using MCMC. These models are applicable in mainly two settings. First, when focus
is on the survival outcome and we wish to account for the effect of an endogenous 
(aka internal) time-dependent covariates measured with error. Second, when focus is on the
longitudinal outcome and we wish to correct for nonrandom dropout.

Some basic features of the package are:

* The user can specify her own density function for the longitudinal responses using argument `densLong` (default is the normal pdf). Among others, this allows to fit joint models with categorical and left-censored longitudinal responses and robust joint models with Student's-t error terms. In addition, using the `df.RE` argument, the user can also change the distribution of the random effects from multivariate normal to a multivariate Student's-t with prespecified degrees of freedom

* The user has now the option to define custom transformation functions for the terms of the longitudinal submodel that enter into the linear predictor of the survival submodel (argument `transFun`). For example, interactions terms, nonlinear terms (polynomials, splines), etc.

* The baseline hazard is estimated using B-splines (penalized (default) or regression). 

* Dynamic predictions
    * function `survfitJM.JMbayes()` computes dynamic survival probabilities.
    * function `predict.JMbayes()` computes dynamic predictions for the longitudinal outcome.
    * function `aucJM()` calculates time-dependent AUCs for joint models, and function `rocJM()` calculates the corresponding time-dependent sensitivities and specifities.
    * function `prederrJM()` calculates prediction errors for joint models.
