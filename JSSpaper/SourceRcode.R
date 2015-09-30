#########################################################################################
# Description: R script replicating the results in the mansuscript entitled 'The R      #
#              Package JMbayes for Fitting Joint Models for Longitudinal and            #
#              Time-to-Event Data using MCMC'                                           #
# Author: Dimitris Rizopoulos                                                           #
# Last update: 2015-07-29                                                               #
#########################################################################################


###############
# Section 4.1 #
###############

library("JMbayes")
library("lattice")
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")


lmeFit.pbc1 <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
                   random = ~ ns(year, 2) | id)
coxFit.pbc1 <- coxph(Surv(years, status2) ~ drug * age, data = pbc2.id, x = TRUE)
jointFit.pbc1 <- jointModelBayes(lmeFit.pbc1, coxFit.pbc1, timeVar = "year", 
                                 n.iter = 30000)


###############
# Section 4.2 #
###############

# a joint model with Student's-t error terms
dLong <- function (y, eta.y, scale, log = FALSE, data) {
    dgt(x = y, mu = eta.y, sigma = scale, df = 4, log = log)
}
jointFit.pbc2 <- jointModelBayes(lmeFit.pbc1, coxFit.pbc1, timeVar = "year", 
                                 densLong = dLong)

# a joint model for a dichotomous longitudinal outcome
pbc2$serBilirD <- as.numeric(pbc2$serBilir > 1.8)
lmeFit.pbc2 <- glmmPQL(serBilirD ~ year, random = ~ year | id, family = binomial, 
                       data = pbc2)
dLongBin <- function (y, eta.y, scale, log = FALSE, data) {
    dbinom(x = y, size = 1, prob = plogis(eta.y), log = log)
}
jointFit.pbc3 <- jointModelBayes(lmeFit.pbc2, coxFit.pbc1, timeVar = "year", 
                                 densLong = dLongBin)

# a joint model for a left censored longitudinal outcome
pbc2$CensInd <- as.numeric(pbc2$serBilir <= 0.8)
pbc2$serBilir2 <- pbc2$serBilir
pbc2$serBilir2[pbc2$serBilir2 <= 0.8] <- 0.8

censdLong <- function (y, eta.y, scale, log = FALSE, data) {
    log.f <- dnorm(x = y, mean = eta.y, sd = scale, log = TRUE)
    log.F <- pnorm(q = y, mean = eta.y, sd = scale, log.p = TRUE)
    ind <- data$CensInd
    log.dens <- (1 - ind) * log.f + ind * log.F
    if (log) log.dens else exp(log.dens)
}
lmeFit.pbc3 <- lme(log(serBilir2) ~ ns(year, 2), data = pbc2,
                   random = ~ ns(year, 2) | id)
jointFit.pbc4 <- jointModelBayes(lmeFit.pbc3, coxFit.pbc1, timeVar = "year",
                                  densLong = censdLong)


###############
# Section 4.3 #
###############

# a joint model with the current value and current slope term 
dForm <- list(fixed = ~ 0 + dns(year, 2), random = ~ 0 + dns(year, 2), 
              indFixed = 2:3, indRandom = 2:3)
jointFit.pbc12 <- update(jointFit.pbc1, param = "td-both", extraForm = dForm)


# a joint model with a cumulative effect
iForm <- list(fixed = ~ 0 + year + ins(year, 2), random = ~ 0 + year + ins(year, 2), 
              indFixed = 1:3, indRandom = 1:3)
jointFit.pbc13 <- update(jointFit.pbc1, param = "td-extra", extraForm = iForm)

wf <- function (u, parms, t.max) {
    num <- dnorm(x = u, sd = parms)
    den <- pnorm(q = c(0, t.max), sd = parms)
    num / (den[2L] - den[1L])
}
jointFit.pbc13w <- update(jointFit.pbc1, estimateWeightFun = TRUE, 
                          weightFun = wf, priorShapes = list(shape1 = dunif),
                          priors = list(priorshape1 = c(0, 10)))


# a joint model with the shared random effects parameterization
jointFit.pbc14 <- update(jointFit.pbc1, param = "shared-RE", n.iter = 50000)


###############
# Section 4.4 #
###############

# use of transformation functions
tf1 <- function (x, data) cbind(x, "^2" = x*x)
tf2 <- function (x, data) cbind(x, "D-penicil" = x * (data$drug == 'D-penicil'))
jointFit.pbc15 <- update(jointFit.pbc12, transFun = list(value = tf1, extra = tf2))


###################################################################################
###################################################################################


###############
# Section 5.1 #
###############

# Dynamic predictions for Patient 2 for the survival and longitudinal outcomes
ND <- pbc2[pbc2$id == 2, ]
sfit.pbc15 <- survfitJM(jointFit.pbc15, newdata = ND)

plot(sfit.pbc15, estimator = "mean", include.y = TRUE,
     conf.int = TRUE, fill.area = TRUE, col.area = "lightgrey")

Ps.pbc15 <- predict(jointFit.pbc15, ND, type = "Subject", 
                    interval = "confidence", return = TRUE)

last.time <- with(Ps.pbc15, year[!is.na(low)][1])

xyplot(pred + low + upp ~ year, data = Ps.pbc15, type = "l", 
       lty = c(1,2,2), col = c(2,1,1), abline = list(v = last.time, lty = 3), 
       xlab = "Time (years)", ylab = "Predicted log(serum bilirubin)")


# Web interface with shiny
runDynPred()


###############
# Section 5.2 #
###############

# Bayesian Model Averaging
Models <- list(jointFit.pbc1, jointFit.pbc12, jointFit.pbc13, 
               jointFit.pbc14, jointFit.pbc15)

log.p.Mk <- log(rep(1/5, 5))
log.p.Dn.Mk <- sapply(Models, logLik, marginal.thetas = TRUE)
log.p.Dj.Mk <- sapply(Models, marglogLik, newdata = ND[1:5, ])

weightsBMA <- log.p.Dj.Mk + log.p.Dn.Mk + log.p.Mk
weightsBMA <- exp(weightsBMA - mean(weightsBMA, na.rm = TRUE))
weightsBMA <- weightsBMA / sum(weightsBMA, na.rm = TRUE)

survPreds <- lapply(Models, survfitJM, newdata = ND[1:5, ])

survPreds.BMA <- bma.combine(JMlis = survPreds, weights = weightsBMA)
survPreds.BMA


###############
# Section 5.3 #
###############

# Discrimination & Calibration
auc.pbc15 <- aucJM(jointFit.pbc15, pbc2, Tstart = 5, Dt = 2)
roc.pbc15 <- rocJM(jointFit.pbc15, pbc2, Tstart = 5, Dt = 2)
dynC.pbc15 <- dynCJM(jointFit.pbc15, pbc2, Dt = 2)
pe.pbc15 <- prederrJM(jointFit.pbc15, pbc2, Tstart = 5, Thoriz = 7)
ipe.pbc15 <- prederrJM(jointFit.pbc15, pbc2, Tstart = 5, Thoriz = 9, interval = TRUE)

# Validation using 10-fold CV
library("parallel")
set.seed(123)
V <- 10
n <- nrow(pbc2.id)
splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
CrossValJM <- function (i) {
    library("JMbayes")
    pbc2$status2 <- as.numeric(pbc2$status != "alive")
    pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")
    trainingData <- pbc2[!pbc2$id %in% i, ]
    trainingData.id <- trainingData[!duplicated(trainingData$id), ]
    testingData <- pbc2[pbc2$id %in% i, ]
    
    lmeFit.pbc1 <- lme(log(serBilir) ~ ns(year, 2), data = trainingData,
                       random = ~ ns(year, 2) | id)
    coxFit.pbc1 <- coxph(Surv(years, status2) ~ drug * age, data = trainingData.id, 
                         x = TRUE)
    
    dForm <- list(fixed = ~ 0 + dns(year, 2), random = ~ 0 + dns(year, 2), 
                  indFixed = 2:3, indRandom = 2:3)
    tf1 <- function (x, data) cbind(x, "^2" = x*x)
    tf2 <- function (x, data) cbind(x, "drugD-penicil" = x * (data$drug == 'D-penicil'))
    jointFit.pbc15 <- jointModelBayes(lmeFit.pbc1, coxFit.pbc1, timeVar = "year",
                                      param = "td-both", extraForm = dForm, 
                                      transFun = list(value = tf1, extra = tf2))
    
    pe <- prederrJM(jointFit.pbc15, newdata = testingData, Tstart = 5, Thoriz = 7)
    auc <- aucJM(jointFit.pbc15, newdata = testingData, Tstart = 5, Thoriz = 7)
    
    list(pe = pe, auc = auc)
}


cl <- makeCluster(5)
res <- parLapply(cl, splits, CrossValJM)
stopCluster(cl)


# cross-validation estimates of the prediction error and the AUC
mean(sapply(res, function (x) x$pe$prederr))
mean(sapply(res, function (x) x$auc$auc))


##############
# Appendix A #
##############

jointFit.pbc1.10knots <- update(jointFit.pbc1, lng.in.kn = 10L)
jointFit.pbc1.20knots <- update(jointFit.pbc1, lng.in.kn = 20L)

jointFit.pbc15.10knots <- update(jointFit.pbc15, lng.in.kn = 10L)
jointFit.pbc15.20knots <- update(jointFit.pbc15, lng.in.kn = 20L)


##############
# Appendix C #
##############

pbc <- pbc2[c("id", "serBilir", "drug", "year", "years",
              "status2", "spiders")]
pbc$start <- pbc$year
splitID <- split(pbc[c("start", "years")], pbc$id)
pbc$stop <- unlist(lapply(splitID,
                          function (d) c(d$start[-1], d$years[1]) ))
pbc$event <- with(pbc, ave(status2, id,
                           FUN = function (x) c(rep(0, length(x)-1), x[1])))
pbc <- pbc[!is.na(pbc$spiders), ]
pbc <- pbc[pbc$start != 0, ]

lmeFit.pbc <- lme(log(serBilir) ~ drug * ns(year, 2),
                  random = ~ ns(year, 2) | id, data = pbc)

tdCox.pbc <- coxph(Surv(start, stop, event) ~ drug * spiders + cluster(id),
                   data = pbc, x = TRUE, model = TRUE)

jointFit.pbc <- jointModelBayes(lmeFit.pbc, tdCox.pbc, timeVar = "year")

summary(jointFit.pbc)

###################################################################################
###################################################################################

save.image(file = "C:/Users/dimitris/Documents/Papers/Paper18/Results/Application/AnalysesPaper.RData")
