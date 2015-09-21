MCMCfit <- function (y, x, param, extraForm, baseHaz, estimateWeightFun, initials, priors, 
                     scales, Funs, Covs, Data, control, df.RE) {
    # extract data longitudinal
    performHC <- control$performHC
    indBetas <- if (performHC) dropAttr(y$indBetas)
    flat_indBetas <- unlist(indBetas)
    y.long <- dropAttr(y$y)
    X <- if (performHC) dropAttr(x$X, flat_indBetas) else dropAttr(x$X)
    Z <- dropAttr(x$Z)
    # extract data survival
    Time <- dropAttr(y$Time)
    event <- dropAttr(y$event)
    TimeL <- dropAttr(y$TimeL)
    W <- dropAttr(x$W)
    notNullW <- !is.null(W)
    Ws <- dropAttr(x$Ws)
    W2 <- dropAttr(x$W2)
    W2s <- dropAttr(x$W2s)
    Xtime <- if (performHC) dropAttr(x$Xtime, flat_indBetas) else dropAttr(x$Xtime)
    XXtime <- dropAttr(x$XXtime)
    Ztime <- dropAttr(x$Ztime)
    iF <- dropAttr(extraForm$indFixed)
    iR <- dropAttr(extraForm$indRandom)
    Xtime.extra <- if (performHC) dropAttr(x$Xtime.extra, flat_indBetas, iF) else dropAttr(x$Xtime.extra)
    Ztime.extra <- dropAttr(x$Ztime.extra)
    Xs <- if (performHC) dropAttr(x$Xs, flat_indBetas) else dropAttr(x$Xs)
    Zs <- dropAttr(x$Zs)
    Xs.extra <- if (performHC) dropAttr(x$Xs.extra, flat_indBetas, iF) else dropAttr(x$Xs.extra)
    Zs.extra <- dropAttr(x$Zs.extra)
    Xu <- if (performHC) dropAttr(x$Xu, flat_indBetas) else dropAttr(x$Xu)
    Zu <- dropAttr(x$Zu)
    # extract indices
    n <- length(Time)
    ns <- nrow(W2s)
    nu <- nrow(Xu)
    ncZ <- ncol(Z)
    nrZ <- nrow(Z)
    nrZtime <- nrow(Ztime)
    nrZs <- nrow(Zs)
    ncZ.extra <- ncol(Ztime.extra)
    nrZtime.extra <- nrow(Ztime.extra)
    nrZs.extra <- nrow(Zs.extra)
    nrZu <- nrow(Zu)
    id <- dropAttr(y$id)
    idFast <- c(id[-length(id)] != id[-1L], TRUE)
    id.GK <- dropAttr(y$id.GK)
    id.GKFast <- c(id.GK[-length(id.GK)] != id.GK[-1L], TRUE)
    idT <- dropAttr(y$idT)
    lag <- y$lag
    LongFormat <- y$LongFormat
    anyLeftCens <- y$anyLeftCens
    typeSurvInfCount <- y$typeSurvInf == "counting"
    w <- rep(dropAttr(x$wk), n)
    P <- dropAttr(x$P)
    st <- c(t(dropAttr(x$st)))
    if (estimateWeightFun) {
        nshapes <- length(initials$shapes)
        seq.nshapes <- seq_len(nshapes)
        weightFun <- Funs$weightFun
        id.GK2 <- dropAttr(y$id.GK2)
        id.GK2Fast <- c(id.GK2[-length(id.GK2)] != id.GK2[-1L], TRUE)
        id.GKu <- rep(id.GK, each = length(x$wk))
        w2 <- rep(dropAttr(x$wk), ns)
        P2 <- dropAttr(x$P2)
        st2 <- c(t(dropAttr(x$st2)))
        max.time <- max(Time)
        u.idGK <- Time[id.GK] - st
        u.idGK2 <- st[id.GK2] - st2
    }
    paramValue <- (param %in% c("td-value", "td-both")) && !estimateWeightFun
    paramExtra <- param %in% c("td-extra", "td-both")
    paramRE <- param %in% c("shared-betasRE", "shared-RE")
    paramSharedRE <- param == "shared-RE"
    baseHazP <- baseHaz == "P-splines"
    paramValueRE <- (paramValue || paramRE)
    estimateAlphas <- paramValueRE || estimateWeightFun
    notestimateWeightFun <- !estimateWeightFun
    rescale_Bs.gammas <- control$rescale_Bs.gammas
    robust_baseHaz <- control$robust_baseHaz
    # extract initial values
    init.betas <- betas <- dropAttr(initials$betas)
    init.tau <- tau <- dropAttr(initials$tau)
    init.b <- b <- dropAttr(initials$b)
    if (param == "shared-betasRE") {
        indBetasRE <- rep(y$indBetasRE, each = n)
    }
    if (performHC) {
        betas1 <- betas[-flat_indBetas]
        betas2 <- betas[flat_indBetas]
        nbetas1 <- length(betas1)
        nbetas2 <- length(betas2)
        mean_b <- function (indBetas) {
            sapply(indBetas, function (i) XXtime[, i, drop = FALSE] %*% betas[i])
        }
        b <- b + mean_b(indBetas)
        indBetasL <- vector("logical", length(betas))
        indBetasL[flat_indBetas] <- TRUE
        XXtimeL <- XXtime[rep(seq_len(n), each = ncZ), ]
        indBetas_ <- lapply(indBetas, function (i) seq_len(nbetas2)[-i])
        row.ind <- rep(seq_len(n * ncZ), rep(sapply(indBetas_, length), n))
        col.ind <- unlist(rep(indBetas_, n))
        XXtimeL[cbind(row.ind, col.ind)] <- 0
        ind <- matrix(seq_len(n * ncZ), ncZ, n)
        XXtimeL2 <- lapply(1:n, function (m) XXtimeL[ind[, m], indBetasL, drop = FALSE])
        XXtimeL3 <- lapply(1:n, function (m) t(XXtimeL[ind[, m], indBetasL, drop = FALSE]))
        nbetas22 <- nbetas2 * nbetas2
        seqn <- seq_len(n)
        diagbetas2 <- diag(nbetas2)
        if (!is.null(iF)) {
            iii <- vector("logical", length(betas))
            iii[iF] <- TRUE
            iF <- which(iii[-flat_indBetas])
        }
    } else {
        betas1 <- betas
        betas2 <- numeric(ncZ)
        nbetas1 <- length(betas1)
        nbetas2 <- length(betas2)
        indBetasL <- vector("logical", length(betas))
    }
    init.invD <- invD <- dropAttr(initials$invD)
    init.gammas <- gammas <- dropAttr(initials$gammas)
    init.Bs.gammas <- Bs.gammas <- dropAttr(initials$Bs.gammas)
    init.tauBs <- tauBs <- dropAttr(initials$tauBs)
    init.deltaBs <- deltaBs <- dropAttr(initials$deltaBs)
    init.alphas <- alphas <- dropAttr(initials$alphas)
    init.Dalphas <- Dalphas <- dropAttr(initials$Dalphas)
    init.shapes <- shapes <- dropAttr(initials$shapes)
    # dimensions of parameters
    nRE <- rep(ncZ, n)
    ngammas <- length(gammas)
    nBs.gammas <- length(Bs.gammas)
    nalphas <- length(alphas)
    nDalphas <- length(Dalphas)
    # extract Funs
    densLong <- Funs$densLong
    hasScale <- Funs$hasScale
    densRE <- Funs$densRE
    transFun.value <- Funs$transFun.value
    transFun.extra <- Funs$transFun.extra
    # Data sets
    data <- Data$data
    data.id <- Data$data.id
    data.s <- Data$data.s
    data.u <- Data$data.u
    # define priors
    priorMean.betas1 <- priors$priorMean.betas[!indBetasL]
    priorTau.betas1 <- priors$priorTau.betas[!indBetasL, !indBetasL]
    log.prior.betas1 <- function (betas1) {
        dmvnorm(betas1, priorMean.betas1, invSigma = priorTau.betas1, log = TRUE)
    }
    priorMean.betas2 <- priors$priorMean.betas[indBetasL]
    priorTau.betas2 <- priors$priorTau.betas[indBetasL, indBetasL]
    log.prior.betas2 <- function (betas2) {
        dmvnorm(betas2, priorMean.betas2, invSigma = priorTau.betas2, log = TRUE)
    }
    priorA.tau <- priors$priorA.tau
    priorB.tau <- priors$priorB.tau
    log.prior.tau <- function (tau) {
        dgamma(tau, priorA.tau, priorB.tau)
    }
    priorR.invD <- priors$priorR.invD
    priorK.invD <- priors$priorK.invD
    log.prior.invD <- function (invD) {
        dwish(invD, priorR.invD, priorK.invD, log = TRUE)
    }
    priorMean.gammas <- priors$priorMean.gammas
    priorTau.gammas <- priors$priorTau.gammas
    log.prior.gammas <- function (gammas) {
        dmvnorm(gammas, priors$priorMean.gammas, invSigma = priorTau.gammas, log = TRUE)
    }
    priorMean.Bs.gammas <- priors$priorMean.Bs.gammas
    priorTau.Bs.gammas <- priors$priorTau.Bs.gammas
    log.prior.Bs.gammas <- function (Bs.gammas) {
        if (baseHazP)
            priorTau.Bs.gammas <- tauBs * priorTau.Bs.gammas
        dmvnorm(Bs.gammas, priorMean.Bs.gammas, invSigma = priorTau.Bs.gammas, log = TRUE)
    }
    priorA.tauBs <- priors$priorA.tauBs
    priorB.tauBs <- priors$priorB.tauBs
    priorA.deltaBs <- priors$priorA.deltaBs
    priorB.deltaBs <- priors$priorB.deltaBs    
    priorMean.alphas <- priors$priorMean.alphas
    priorTau.alphas <- priors$priorTau.alphas
    log.prior.alphas <- function (alphas) {
        dmvnorm(alphas, priorMean.alphas, invSigma = priorTau.alphas, log = TRUE)
    }
    priorMean.Dalphas <- priors$priorMean.Dalphas
    priorTau.Dalphas <- priors$priorTau.Dalphas
    log.prior.Dalphas <- function (Dalphas) {
        dmvnorm(Dalphas, priorMean.Dalphas, invSigma = priorTau.Dalphas, log = TRUE)
    }
    priorshape1Fun <- control$priorShapes$shape1
    priorshape1.1 <- priors$priorshape1[1L]
    priorshape1.2 <- priors$priorshape1[2L]
    log.prior.shape1 <- function (shape1) {
        priorshape1Fun(shape1, priorshape1.1, priorshape1.2, log = TRUE)
    }
    priorshape2Fun <- control$priorShapes$shape2
    priorshape2.1 <- priors$priorshape2[1L]
    priorshape2.2 <- priors$priorshape2[2L]
    log.prior.shape2 <- function (shape2) {
        priorshape2Fun(shape2, priorshape2.1, priorshape2.2, log = TRUE)
    }
    priorshape3Fun <- control$priorShapes$shape3
    priorshape3.1 <- priors$priorshape3[1L]
    priorshape3.2 <- priors$priorshape3[2L]
    log.prior.shape3 <- function (shape3) {
        priorshape3Fun(shape3, priorshape3.1, priorshape3.2, log = TRUE)
    }
    # define posteriors
    logPost.betas1 <- function (betas1) {
        Xbetas <- drop(X %*% betas1)
        eta.y <- Xbetas + Zb
        log.pyb <- fastSumID2(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), idFast)
        log.prior <- log.prior.betas1(betas1)
        if (!paramRE) {
            Mtime <- numeric(n)
            Ms <- numeric(ns)
            if (paramValue) {
                Xtimebetas <- drop(Xtime %*% betas1)
                Xsbetas <- drop(Xs %*% betas1)
                vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
                vls <- transFun.value(Xsbetas + Zsb, data.s)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Xtime.extrabetas <- drop(Xtime.extra %*% betas1[iF])
                Xs.extrabetas <- drop(Xs.extra %*% betas1[iF])
                ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
                exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            if (estimateWeightFun) {
                Xsbetas <- drop(Xs %*% betas1)
                Xubetas <- drop(Xu %*% betas1)
                vl <- transFun.value(P * fastSumID2(wFun * (Xsbetas + Zsb), id.GKFast), data.id)
                vls <- transFun.value(P2 * fastSumID2(wFun2 * (Xubetas + Zub), id.GK2Fast), data.s)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (notNullW) {
                if (LongFormat) {
                    log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
                } else {
                    Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                    log.Surv <- expWgammas * Int
                }
            } else {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
            }
            log.ptb <- event * Mtime - log.Surv
            list(log.post = sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior,
                 Xbetas = Xbetas,
                 Xtimebetas = if (paramValue) Xtimebetas,
                 Xsbetas = if (estimateAlphas) Xsbetas,
                 Xtime.extrabetas = if (paramExtra) Xtime.extrabetas,
                 Xs.extrabetas = if (paramExtra) Xs.extrabetas,
                 Xubetas = if (estimateWeightFun) Xubetas,
                 vl = if (estimateAlphas) vl, vls = if (estimateAlphas) vls,
                 ex = if (paramExtra) ex, exs = if (paramExtra) exs, 
                 Ms = Ms, Mtime = Mtime, log.Surv = log.Surv, Int = Int)
        } else {
            if (paramSharedRE) {
                list(log.post = sum(log.pyb, na.rm = TRUE) + log.prior, Xbetas = Xbetas)
            }
        }
    }
    logPost.betas1Fast <- function () {
        log.pyb <- fastSumID2(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), idFast)
        log.prior <- log.prior.betas1(betas1)
        if (!paramRE) {
            Mtime <- numeric(n)
            if (paramValue) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            if (paramExtra) {
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            }
            if (estimateWeightFun) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            log.ptb <- event * Mtime - log.Surv
            sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior
        } else {
            if (paramSharedRE) {
                sum(log.pyb, na.rm = TRUE) + log.prior
            }
        }
    }
    logPost.tau <- function (tau) {
        log.pyb <- densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data)
        log.prior <- log.prior.tau(tau)
        sum(log.pyb, na.rm = TRUE) + log.prior
    }
    logPost.RE <- function (b) {
        Zb <- .rowSums(Z * b[id, , drop = FALSE], nrZ, ncZ)
        eta.y <- if (nbetas1) Xbetas + Zb else Zb
        log.pyb <- fastSumID2(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), idFast)
        mu_b <- if (performHC) mean_b(indBetas) else betas2
        log.prior <- densRE(b, mu = mu_b, invD = invD, log = TRUE)
        if (!paramRE) {
            Mtime <- numeric(n)
            Ms <- numeric(ns)
            if (paramValue) {
                Zsb <- .rowSums(Zs * b[id.GK, , drop = FALSE], nrZs, ncZ)
                Ztimeb <- .rowSums(Ztime * b, nrZtime, ncZ)
                if (nbetas1) {
                    vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
                    vls <- transFun.value(Xsbetas + Zsb, data.s)                    
                } else {
                    vl <- transFun.value(Ztimeb, data.id)
                    vls <- transFun.value(Zsb, data.s)                    
                }
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Ztime.extrab <- .rowSums(Ztime.extra * b[, iR, drop = FALSE], nrZtime.extra, ncZ.extra)
                Zs.extrab <- .rowSums(Zs.extra * b[id.GK, iR, drop = FALSE], nrZs.extra, ncZ.extra)
                if (nbetas1) {
                    ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
                    exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)                    
                } else {
                    ex <- transFun.extra(Ztime.extrab, data.id)
                    exs <- transFun.extra(Zs.extrab, data.s)
                }
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            if (estimateWeightFun) {
                Zsb <- .rowSums(Zs * b[id.GK, , drop = FALSE], nrZs, ncZ)
                Zub <- .rowSums(Zu * b[id.GKu, , drop = FALSE], nrZu, ncZ)
                if (nbetas1) {
                    vl <- transFun.value(P * fastSumID2(wFun * (Xsbetas + Zsb), id.GKFast), data.id)
                    vls <- transFun.value(P2 * fastSumID2(wFun2 * (Xubetas + Zub), id.GK2Fast), data.s)
                } else {
                    vl <- transFun.value(P * fastSumID2(wFun * Zsb, id.GKFast), data.id)
                    vls <- transFun.value(P2 * fastSumID2(wFun2 * Zub, id.GK2Fast), data.s)
                }
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }            
            if (notNullW) {
                if (LongFormat) {
                    log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
                } else {
                    Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                    log.Surv <- expWgammas * Int
                }
            } else {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
            }
            log.ptb <- event * Mtime - log.Surv
            list(log.post = log.pyb + log.ptb + log.prior,
                 Zb = Zb, Int = Int, Ms = Ms, eta.y = eta.y, log.Surv = log.Surv,
                 Ztimeb = if (paramValue) Ztimeb,
                 Zsb = if (estimateAlphas) Zsb,
                 Ztime.extrab = if (paramExtra) Ztime.extrab,
                 Zs.extrabetas = if (paramExtra) Zs.extrab,
                 Zub = if (estimateWeightFun) Zub,
                 vl = if (estimateAlphas) vl,
                 vls = if (estimateAlphas) vls,
                 ex = if (paramExtra) ex,
                 exs = if (paramExtra) exs)
        } else {
            Mtime <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
            log.Surv <- exp(Mtime) * Int
            if (notNullW && !LongFormat)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            list(log.post = log.pyb + log.ptb + log.prior,
                 Zb = Zb, eta.y = eta.y, log.Surv = log.Surv)
        }
    }
    logPost.REFast <- function () {
        eta.y <- if (nbetas1) Xbetas + Zb else Zb
        log.pyb <- fastSumID2(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), idFast)
        mu_b <- if (performHC) mean_b(indBetas) else betas2
        log.prior <- densRE(b, mu = mu_b, invD = invD, log = TRUE)
        if (!paramRE) {
            Mtime <- numeric(n)
            if (paramValue) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            if (paramExtra) {
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            }
            if (estimateWeightFun) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            log.ptb <- event * Mtime - log.Surv
        } else {
            Mtime <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
            log.ptb <- event * Mtime - log.Surv
        }
        log.pyb + log.ptb + log.prior
    }
    logPost.invD <- function (invD) {
        log.pb <- densRE(b, invD = invD, log = TRUE, prop = FALSE)
        log.prior <- log.prior.invD(invD)
        sum(log.pb, na.rm = TRUE) + log.prior
    }
    logPost.gammas <- function (gammas){
        Wgammas <- drop(W %*% gammas)
        log.Surv <- if (LongFormat) {
            Wsgammas <- drop(Ws %*% gammas)
            log.integrand <- if (paramRE) log.h0s + Wsgammas else log.h0s + Wsgammas + Ms
            Int <- P * fastSumID2(w * exp(log.integrand), id.GKFast)
            Int
        } else {
            exp(Wgammas) * Int
        }
        if (paramRE) {
            Mtime <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
            log.Surv <- exp(Mtime) * log.Surv
        }
        log.ptb <- event * Wgammas - log.Surv
        log.prior <- log.prior.gammas(gammas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, expWgammas = exp(Wgammas),
             log.Surv = log.Surv, Int = if (LongFormat) Int, 
             Wsgammas = if (LongFormat) Wsgammas)
    }
    logPost.Bs.gammas <- function (Bs.gammas) {
        if (rescale_Bs.gammas)
            Bs.gammas <- drop(tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas)
        W2sBs.gammas <- drop(W2s %*% Bs.gammas)
        log.integrand <- if (paramRE) W2sBs.gammas else W2sBs.gammas + Ms
        if (notNullW && LongFormat) {
            log.integrand <- log.integrand + Wsgammas
        }
        log.Surv <- Int <- P * fastSumID2(w * exp(log.integrand), id.GKFast) 
        if (paramRE) {
            Mtime <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
            log.Surv <- Int * exp(Mtime)
        }
        if (notNullW && !LongFormat) {
            log.Surv <- expWgammas * log.Surv
        }
        log.ptb <- event * drop(W2 %*% Bs.gammas) - log.Surv
        log.prior <- log.prior.Bs.gammas(Bs.gammas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.h0s = W2sBs.gammas,
             log.Surv = log.Surv, Int = Int)
    }
    logPost.Bs.gammasFast <- function () {
        if (rescale_Bs.gammas)
            Bs.gammas <- drop(tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas)
        log.ptb <- event * drop(W2 %*% Bs.gammas) - log.Surv
        log.prior <- log.prior.Bs.gammas(Bs.gammas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    ArankDiff <- priorA.tauBs + 0.5 * qr(priorTau.Bs.gammas)$rank
    AdeltaBs <- priorA.deltaBs + priorA.tauBs
    logPost.alphas <- function (alphas) {
        if (!paramRE) {
            Ms <- numeric(ns)
            if (estimateAlphas) {
                Mtime.alphas <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            if (notNullW) {
                if (LongFormat) {
                    log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
                } else {
                    Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                    log.Surv <- expWgammas * Int
                }
            } else {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
            }            
            log.ptb <- event * Mtime.alphas - log.Surv
            log.prior <- log.prior.alphas(alphas)
            list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv, 
                 Ms = Ms, Int = Int)
        } else {
            Mtime <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
            log.Surv <- exp(Mtime) * Int
            if (notNullW && !LongFormat)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            log.prior <- log.prior.alphas(alphas)
            list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv)
        }
    }
    logPost.alphasFast <- function () {
        if (!paramRE) {
            Mtime.alphas <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        } else {
            Mtime.alphas <- if (paramSharedRE) {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                drop((b - mu_b) %*% alphas) 
            } else {
                mu_b <- if (performHC) mean_b(indBetas) else betas2
                betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                if (nbetas1) betas[!indBetasL] <- betas1
                if (performHC) betas[indBetasL] <- betas2
                bb_ <- (b - mu_b) + betas[indBetasRE]
                drop(bb_ %*% alphas)
            }
        }
        log.ptb <- event * Mtime.alphas - log.Surv
        log.prior <- log.prior.alphas(alphas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    logPost.Dalphas <- function (Dalphas) {
        Ms <- numeric(ns)
        if (estimateAlphas) {
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        if (paramExtra) {
            Mtime.Dalphas <- if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }        
        if (notNullW) {
            if (LongFormat) {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
            } else {
                Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                log.Surv <- expWgammas * Int
            }
        } else {
            log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
        }
        log.ptb <- event * Mtime.Dalphas - log.Surv
        log.prior <- log.prior.Dalphas(Dalphas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv,
             Ms = Ms, Int = Int)
    }
    logPost.DalphasFast <- function () {
        Mtime.Dalphas <- if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
        log.ptb <- event * Mtime.Dalphas - log.Surv
        log.prior <- log.prior.Dalphas(Dalphas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    logPost.shape <- function (shape, which) {
        shapes[which] <- shape
        Ms <- numeric(ns)
        ###
        wFun <- w * weightFun(u.idGK, shapes, max.time)
        wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
        if (nbetas1) {
            vl <- transFun.value(P * fastSumID2(wFun * XsbetasZsb, id.GKFast), data.id)
            vls <- transFun.value(P2 * fastSumID2(wFun2 * XubetasZub, id.GK2Fast), data.s)
        } else {
            vl <- transFun.value(P * fastSumID2(wFun * Zsb, id.GKFast), data.id)
            vls <- transFun.value(P2 * fastSumID2(wFun2 * Zub, id.GK2Fast), data.s)
        }
        Mtime <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        ###
        if (paramExtra) {
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }
        ###
        if (notNullW) {
            if (LongFormat) {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
            } else {
                Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                log.Surv <- expWgammas * Int
            }
        } else {
            log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
        }
        log.ptb <- event * Mtime - log.Surv
        log.prior <- switch(which, "1" = log.prior.shape1(shape), 
                            "2" = log.prior.shape2(shape), "3" = log.prior.shape3(shape))
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv,
             Ms = Ms, Int = Int, wFun = wFun, wFun2 = wFun2, vl = vl, vls = vls)
    }
    # define proposals
    # betas
    if (nbetas1) {
        propCov.betas1 <- eigen(Covs$betas[!indBetasL, !indBetasL], symmetric = TRUE)
        scale.betas1 <- if (!is.null(ss <- scales[["betas"]])) ss else 5.66/nbetas1
        r.betas1 <- function (n) {
            propCov.betas1$values <- propCov.betas1$values * scale.betas1
            rmvnorm(n, mu = NULL, Sigma = propCov.betas1)
        }        
    }
    # b
    propCov.RE <- lapply(Covs$b, eigen, symmetric = TRUE)
    scale.RE <- if (!is.null(ss <- scales[["b"]])) ss else 5.66/nRE
    r.RE <- function (N) {
        out <- array(0, c(dim(b), N))
        for (i in 1:n) {
            propCov.RE[[i]][["values"]] <- propCov.RE[[i]][["values"]] * scale.RE[i]
            out[i, , ] <- rmvnorm(N, mu = NULL, Sigma = propCov.RE[[i]])
        }
        out
    }
    # invD
    isNulldf.RE <- is.null(df.RE)
    diagB <- diag(1, ncZ)
    if (isNulldf.RE) {
        K.invDn <- priorK.invD + n
        r.invD <- function (N) {
            drop(rWishart(N, K.invDn, R.Dbtb))
        }
    } else {
        K.invDn <- (df.RE / (df.RE - 2)) * (priorK.invD + n)
        r.invD <- function (N) {
            drop(rWishart(N, K.invDn, invD / K.invDn))
        }
    }
    # gammas
    if (notNullW) {
        propCov.gammas <- eigen(Covs$gammas, symmetric = TRUE)
        scale.gammas <- if (!is.null(ss <- scales$gammas)) ss else 5.66/ngammas
        r.gammas <- function (N) {
            propCov.gammas$values <- propCov.gammas$values * scale.gammas
            rmvnorm(N, mu = NULL, Sigma = propCov.gammas)
        }
    }
    # Bs.gammas
    if (rescale_Bs.gammas) {
        tchol_CovBs.gammas <- t(chol(Covs$Bs.gammas))
        Bs.gammas <- rep(0, nBs.gammas)
        scale.Bs.gammas <- if (!is.null(ss <- scales$Bs.gammas)) ss else 5.66/nBs.gammas
        r.Bs.gammas <- function (N) {
            matrix(rnorm(N * nBs.gammas, sd = scale.Bs.gammas), N, nBs.gammas)
        }
    } else {
        propCov.Bs.gammas <- eigen(Covs$Bs.gammas, symmetric = TRUE)
        scale.Bs.gammas <- if (!is.null(ss <- scales$Bs.gammas)) ss else 5.66/nBs.gammas
        r.Bs.gammas <- function (N) {
            propCov.Bs.gammas$values <- propCov.Bs.gammas$values * scale.Bs.gammas
            rmvnorm(N, mu = NULL, Sigma = propCov.Bs.gammas)
        }
    }
    # alphas
    if (estimateAlphas) {
        propCov.alphas <- eigen(Covs$alphas, symmetric = TRUE)
        scale.alphas <- if (!is.null(ss <- scales$alphas)) ss else 5.66/nalphas
        r.alphas <- function (N) {
            propCov.alphas$values <- propCov.alphas$values * scale.alphas
            rmvnorm(N, mu = NULL, Sigma = propCov.alphas)
        }
    }
    # Dalphas
    if (paramExtra) {
        propCov.Dalphas <- eigen(Covs$Dalphas, symmetric = TRUE)
        scale.Dalphas <- if (!is.null(ss <- scales$Dalphas)) ss else 5.66/nDalphas
        r.Dalphas <- function (N) {
            propCov.Dalphas$values <- propCov.Dalphas$values * scale.Dalphas
            rmvnorm(N, mu = NULL, Sigma = propCov.Dalphas)
        }
    }
    # number of iterations
    n.adapt <- control$n.adapt
    n.burnin <- control$n.burnin
    totalIter <- control$n.iter + n.adapt + n.burnin
    n.thin <- control$n.thin
    n.batch <- control$n.batch
    # objects to keep results
    resInd <- seq(n.adapt + n.burnin + 1L, totalIter, by = n.thin)
    n.out <- length(resInd)
    if (nbetas1)
        res.betas1 <- matrix(0, n.out, nbetas1)
    if (performHC)
        res.betas2 <- matrix(0, n.out, nbetas2)
    if (hasScale)
        res.tau <- matrix(0, n.out, 1)
    res.b <- array(0, c(dim(b), n.out))
    if (performHC)
        res.mean_b <- array(0, c(dim(b), n.out))
    res.invD <- matrix(0, n.out, length(invD))
    res.Bs.gammas <- matrix(0, n.out, length(Bs.gammas))
    if (baseHazP)
        res.tauBs <- matrix(0, n.out, 1)
    if (notNullW)
        res.gammas <- matrix(0, n.out, length(gammas))
    if (estimateAlphas)
        res.alphas <- matrix(0, n.out, length(alphas))
    if (paramExtra)
        res.Dalphas <- matrix(0, n.out, length(Dalphas))
    if (estimateWeightFun)
        res.shapes <- matrix(0, n.out, length(shapes))
    res.logLik <- matrix(0, n.out, n)
    # acceptance rates
    ar.betas <- ar.invD <- ar.gammas <- ar.Bs.gammas <- ar.alphas <- ar.Dalphas <- numeric(totalIter)
    ar.b <- matrix(0, totalIter, n)
    # initiate all components at the starting values
    Zb <- rowSums(Z * b[id, , drop = FALSE])
    eta.y <- if (nbetas1) {
        Xbetas <- drop(X %*% betas1)
        Xbetas + Zb 
    } else Zb
    Mtime <- numeric(n)
    Ms <- numeric(ns)
    if (paramValue) {
        Ztimeb <- rowSums(Ztime * b)
        Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
        if (nbetas1) {
            Xtimebetas <- drop(Xtime %*% betas1)
            Xsbetas <- drop(Xs %*% betas1)
            vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
            vls <- transFun.value(Xsbetas + Zsb, data.s)
        } else {
            vl <- transFun.value(Ztimeb, data.id)
            vls <- transFun.value(Zsb, data.s)
        }
        is.matrix.vl <- is.matrix(vl); is.matrix.vls <- is.matrix(vls)
        Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
    }
    if (paramExtra) {
        Ztime.extrab <- rowSums(Ztime.extra * b[, iR, drop = FALSE])
        Zs.extrab <- rowSums(Zs.extra * b[id.GK, iR, drop = FALSE])
        if (nbetas1) {
            Xtime.extrabetas <- drop(Xtime.extra %*% betas1[iF])
            Xs.extrabetas <- drop(Xs.extra %*% betas1[iF])
            ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
            exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)
        } else {
            ex <- transFun.extra(Ztime.extrab, data.id)
            exs <- transFun.extra(Zs.extrab, data.s)
        }
        is.matrix.ex <- is.matrix(ex); is.matrix.exs <- is.matrix(exs)
        Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
        Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
    }
    if (estimateWeightFun) {
        wFun <- w * weightFun(u.idGK, shapes, max.time)
        Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
        wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
        Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
        if (nbetas1) {
            Xsbetas <- drop(Xs %*% betas1)
            Xubetas <- drop(Xu %*% betas1)
            vl <- transFun.value(P * fastSumID2(wFun * (Xsbetas + Zsb), id.GKFast), data.id)
            vls <- transFun.value(P2 * fastSumID2(wFun2 * (Xubetas + Zub), id.GK2Fast), data.s)
        } else {
            vl <- transFun.value(P * fastSumID2(wFun * Zsb, id.GKFast), data.id)
            vls <- transFun.value(P2 * fastSumID2(wFun2 * Zub, id.GK2Fast), data.s)
        }
        is.matrix.vl <- is.matrix(vl); is.matrix.vls <- is.matrix(vls)
        Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
    }
    log.h0s <- if (rescale_Bs.gammas) {
        drop(W2s %*% (tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas))
    } else drop(W2s %*% Bs.gammas)
    log.integrand <- if (paramRE) log.h0s else log.h0s + Ms
    if (notNullW && LongFormat) {
        Wsgammas <- drop(Ws %*% gammas)
        log.integrand <- log.integrand + Wsgammas
    }
    log.Surv <- Int <- P * fastSumID2(w * exp(log.integrand), id.GKFast) 
    if (paramRE) {
        Mtime <- if (paramSharedRE) {
            mu_b <- if (performHC) mean_b(indBetas) else betas2
            drop((b - mu_b) %*% alphas) 
        } else {
            mu_b <- if (performHC) mean_b(indBetas) else betas2
            betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
            if (nbetas1) betas[!indBetasL] <- betas1
            if (performHC) betas[indBetasL] <- betas2
            bb_ <- (b - mu_b) + betas[indBetasRE]
            drop(bb_ %*% alphas)
        }
        log.Surv <- Int * exp(Mtime)
    }
    if (notNullW && !LongFormat) {
        expWgammas <- exp(drop(W %*% gammas))
        log.Surv <- expWgammas * log.Surv
    }
    #########################################################################################################
    # run the MCMC
    set.seed(control$seed)
    batch <- 1L
    jj <- 0L
    batchStart <- round((n.adapt + n.burnin) / n.batch)
    if (control$verbose) {
        cat("\n MCMC iterations:\n\n")
        pb <- txtProgressBar(0, totalIter, style = 3, char = "+", width = 50)
    }
    time <- system.time(for (i in seq_len(totalIter)) {#
        if (i == 1L || !i %% n.batch) {
            if (i > 1 && i <= n.adapt) {
                ss <- seq(1L + n.batch * (batch - 1L), n.batch * batch)
                if (nbetas1 && is.null(scales[["betas"]]))
                    scale.betas1 <- scbetasF <- adjustScaleRW(scale.betas1, mean(ar.betas[ss]), nbetas1)
                if (is.null(scales[["b"]]))
                    scale.RE <- mapply(adjustScaleRW, scale = scale.RE, acceptRate = colMeans(ar.b[ss, ]), d = nRE)
                if (notNullW && is.null(scales$gammas)) {
                    scale.gammas <- scgammasF <- adjustScaleRW(scale.gammas, mean(ar.gammas[ss]), ngammas)
                }
                if (is.null(scales$Bs.gammas))
                    scale.Bs.gammas <- scBs.gammasF <- adjustScaleRW(scale.Bs.gammas, mean(ar.Bs.gammas[ss]), nBs.gammas)
                if (estimateAlphas && is.null(scales$alphas)) {
                    scale.alphas <- scalphasF <- adjustScaleRW(scale.alphas, mean(ar.alphas[ss]), nalphas)
                }
                if (paramExtra && is.null(scales$Dalphas)) {
                    scale.Dalphas <- adjustScaleRW(scale.Dalphas, mean(ar.Dalphas[ss]), nDalphas)
                }
                if (!isNulldf.RE) {
                    K.invDn <- adjustKRW(K.invDn, mean(ar.invD[ss]), ncZ)
                }
            }
            if (i > 1)
                batch <- batch + 1L
            if (nbetas1)
                new.betas1 <- r.betas1(n.batch)
            new.b <- r.RE(n.batch)
            if (notNullW)
                new.gammas <- r.gammas(n.batch)
            new.Bs.gammas <- r.Bs.gammas(n.batch)
            if (estimateAlphas)
                new.alphas <- r.alphas(n.batch)
            if (paramExtra)
                new.Dalphas <- r.Dalphas(n.batch)
        }
        # batch index
        ii <- i - n.batch * if (!i %% n.batch) i %/% n.batch - 1L else i %/% n.batch
        # update betas1
        if (nbetas1) {
            lP.old.betas <- logPost.betas1Fast()
            new.betas1[ii, ] <- new.betas1[ii, ] + betas1
            lP.betas <- logPost.betas1(new.betas1[ii, ])
            lP.new.betas <- lP.betas$log.post
            lRatio.betas <- lP.new.betas - lP.old.betas
            if (is.finite(lRatio.betas) && (lRatio.betas >= 0 || runif(1L) < exp(lRatio.betas))) {
                ar.betas[i] <- 1
                betas1 <- new.betas1[ii, ]
                Xbetas <- lP.betas$Xbetas
                eta.y <- Xbetas + Zb
                if (!paramRE) {
                    Xtimebetas <- lP.betas$Xtimebetas
                    Xsbetas <- lP.betas$Xsbetas
                    Xtime.extrabetas <- lP.betas$Xtime.extrabetas
                    Xs.extrabetas <- lP.betas$Xs.extrabetas
                    Xubetas <- lP.betas$Xubetas
                    vl <- lP.betas$vl; vls <- lP.betas$vls
                    ex <- lP.betas$ex; exs <- lP.betas$exs
                    Mtime <- lP.betas$Mtime
                    Ms <- lP.betas$Ms
                    log.Surv <- lP.betas$log.Surv
                    Int <- lP.betas$Int
                }
                if (param == "shared-betasRE")
                    log.Surv <- lP.betas$log.Surv
            }
        }
        # update tau
        if (hasScale) {
            tau <- slice.tau(logPost.tau, tau, step = 0.5)
        }
        # update RE
        lP.old.b <- logPost.REFast() #logPost.RE(b)[[1]] 
        new.b[, , ii] <- new.b[, , ii] + b
        lP.RE <- logPost.RE(as.matrix(new.b[, , ii]))
        lP.new.b <- lP.RE$log.post
        lRatio.b <- lP.new.b - lP.old.b
        indRE <- runif(n) < pmin(exp(lRatio.b), 1)
        if (anyNA(indRE))
            indRE[is.na(indRE)] <- FALSE
        indRE.GK <- indRE[id.GK]
        indRE.id <- indRE[id]
        ar.b[i, indRE] <- 1
        b[indRE, ] <- new.b[indRE, , ii]
        Zb[indRE.id] <- lP.RE$Zb[indRE.id]
        log.Surv[indRE] <- lP.RE$log.Surv[indRE]
        eta.y <- if (nbetas1) Xbetas + Zb else Zb
        if (!paramRE) {
            Ms[indRE.GK] <- lP.RE$Ms[indRE.GK]
            if (notNullW) {
                if (LongFormat) {
                    log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
                } else {
                    Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                    log.Surv <- expWgammas * Int
                }
            } else {
                log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
            }
            if (estimateAlphas) {
                if (notestimateWeightFun)
                    Ztimeb[indRE] <- lP.RE$Ztimeb[indRE]
                Zsb[indRE.GK] <- lP.RE$Zsb[indRE.GK]
                if (is.matrix.vl) {
                    vl[indRE, ] <- lP.RE$vl[indRE, ]
                    vls[indRE.GK, ] <- lP.RE$vls[indRE.GK, ]
                } else {
                    vl[indRE] <- lP.RE$vl[indRE]
                    vls[indRE.GK] <- lP.RE$vls[indRE.GK]
                }
            }
            if (paramExtra) {
                Ztime.extrab[indRE] <- lP.RE$Ztime.extrab[indRE]
                Zs.extrab[indRE.GK] <- lP.RE$Zs.extrab[indRE.GK]
                if (is.matrix.ex) {
                    ex[indRE, ] <- lP.RE$ex[indRE, ]
                    exs[indRE.GK, ] <- lP.RE$exs[indRE.GK, ]
                } else {
                    ex[indRE] <- lP.RE$ex[indRE]
                    exs[indRE.GK] <- lP.RE$exs[indRE.GK]
                }
            }
            if (estimateWeightFun) {
                indRE.GKu <- indRE[id.GKu]
                Zub[indRE.GKu] <- lP.RE$Zub[indRE.GKu]
            }
        }
        # update betas2
        if (performHC) {
            mat <- numeric(nbetas22); dim(mat) <- c(nbetas2, nbetas2)
            mu <- numeric(nbetas2)
            for (m in seqn) {
                xxt_invD <- XXtimeL3[[m]] %*% invD
                mat <- mat + xxt_invD %*% XXtimeL2[[m]]
                mu <- mu + xxt_invD %*% b[m, ]
            }
            invV_betas2 <- mat + priorTau.betas2
            V_betas2 <- solve.default(invV_betas2, diagbetas2)
            mu_betas2 <- V_betas2 %*% mu
            betas2 <- rmvnorm(1, mu_betas2, V_betas2)
        }
        # update invD
        if (isNulldf.RE) {
            bb <- if (performHC) crossprod(b - mean_b(indBetas)) else crossprod(b)
            R.Dbtb <- solve.default(priorR.invD + bb, diagB)
            new.invD <- r.invD(1)
            ar.invD[i] <- 1
            invD <- new.invD
        } else {
            new.invD <- r.invD(1)
            lP.old.invD <- logPost.invD(invD)
            lP.new.invD <- logPost.invD(new.invD)
            lRatio.invD <- lP.new.invD + dwish(invD, invD / K.invDn, K.invDn, TRUE) -
                lP.old.invD - dwish(new.invD, invD / K.invDn, K.invDn, TRUE)
            if (lRatio.invD >= 0 || runif(1L) < exp(lRatio.invD)) {
                ar.invD[i] <- 1
                invD <- new.invD
            }
        }
        # update gammas
        if (notNullW) {
            lP.old.gammas <- logPost.gammas(gammas)$log.post
            new.gammas[ii, ] <- new.gammas[ii, ] + gammas
            lP.gammas <- logPost.gammas(new.gammas[ii, ])
            lP.new.gammas <- lP.gammas$log.post
            lRatio.gammas <- lP.new.gammas - lP.old.gammas
            if (is.finite(lRatio.gammas) && (lRatio.gammas >= 0 || runif(1L) < exp(lRatio.gammas))) {
                ar.gammas[i] <- 1
                gammas <- new.gammas[ii, ]
                expWgammas <- lP.gammas$expWgammas
                log.Surv <- lP.gammas$log.Surv
                if (LongFormat) {
                    Wsgammas <- lP.gammas$Wsgammas
                    Int <- lP.gammas$Int
                }
            }
        }
        # update Bs.gammas
        lP.old.Bs.gammas <- logPost.Bs.gammasFast()
        new.Bs.gammas[ii, ] <- new.Bs.gammas[ii, ] + Bs.gammas
        lP.Bs.gammas <- logPost.Bs.gammas(new.Bs.gammas[ii, ])
        lP.new.Bs.gammas <- lP.Bs.gammas$log.post
        lRatio.Bs.gammas <- lP.new.Bs.gammas - lP.old.Bs.gammas
        if (is.finite(lRatio.Bs.gammas) && (lRatio.Bs.gammas >= 0 || runif(1L) < exp(lRatio.Bs.gammas))) {
            ar.Bs.gammas[i] <- 1
            Bs.gammas <- new.Bs.gammas[ii, ]
            log.h0s <- lP.Bs.gammas$log.h0s
            log.Surv <- lP.Bs.gammas$log.Surv
            Int <- lP.Bs.gammas$Int
        }
        # update tauBs
        if (baseHazP) {
            Bs.gammas_s <- if (rescale_Bs.gammas) {
                drop(tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas)
            } else Bs.gammas
            if (robust_baseHaz) {
                BB <- deltaBs * priorB.tauBs + 
                    0.5 * drop(crossprod(Bs.gammas_s, priorTau.Bs.gammas %*% Bs.gammas_s))
                tauBs <- rgamma(1L, ArankDiff, BB)
                deltaBs <- rgamma(1L, AdeltaBs, priorB.deltaBs + priorB.tauBs * tauBs)                
            }
            else {
                BB <- priorB.tauBs + 
                    0.5 * drop(crossprod(Bs.gammas_s, priorTau.Bs.gammas %*% Bs.gammas_s))
                tauBs <- rgamma(1L, ArankDiff, BB)                
            }
        }
        # update alphas
        if (estimateAlphas) {
            lP.old.alphas <- logPost.alphasFast()
            new.alphas[ii, ] <- new.alphas[ii, ] + alphas
            lP.alphas <- logPost.alphas(new.alphas[ii, ])
            lP.new.alphas <- lP.alphas$log.post
            lRatio.alphas <- lP.new.alphas - lP.old.alphas
            if (is.finite(lRatio.alphas) && (lRatio.alphas >= 0 || runif(1L) < exp(lRatio.alphas))) {
                ar.alphas[i] <- 1
                alphas <- new.alphas[ii, ]
                log.Surv <- lP.alphas$log.Surv
                if (!paramRE) {
                    Ms <- lP.alphas$Ms
                    Int <- lP.alphas$Int
                }
            }
        }
        # update Dalphas
        if (paramExtra) {
            lP.old.Dalphas <- logPost.DalphasFast()
            new.Dalphas[ii, ] <- new.Dalphas[ii, ] + Dalphas
            lP.Dalphas <- logPost.Dalphas(new.Dalphas[ii, ])
            lP.new.Dalphas <- lP.Dalphas$log.post
            lRatio.Dalphas <- lP.new.Dalphas - lP.old.Dalphas
            if (is.finite(lRatio.Dalphas) && (lRatio.Dalphas >= 0 || runif(1L) < exp(lRatio.Dalphas))) {
                ar.Dalphas[i] <- 1
                Dalphas <- new.Dalphas[ii, ]
                Ms <- lP.Dalphas$Ms
                log.Surv <- lP.Dalphas$log.Surv
                Int <- lP.Dalphas$Int
            }
        }
        # update shapes
        if (estimateWeightFun) {
            if (control$verbose2)
                cat("\ni =", i, "\tshapes =", round(shapes, 3L), 
                    if (paramExtra) "\tDalphas =", if (paramExtra) round(Dalphas, 3L), 
                    "\talphas =", round(alphas, 3L), "\tbetas = ", round(betas, 3L))
            if (nbetas1) {
                XsbetasZsb <- Xsbetas + Zsb
                XubetasZub <- Xubetas + Zub
            }
            for (shp in seq.nshapes) {
                ss <- 1
                slice.shp <- slice.shape(logPost.shape, shapes, step = ss, which = shp)
                while (slice.shp$fail) {
                    ss <- ss/10
                    if (ss < 1e-03)
                        break
                    slice.shp <- slice.shape(logPost.shape, shapes, step = ss, which = shp)
                }
                if (!slice.shp$fail) {
                    shapes[shp] <- slice.shp$new.shape
                    if (shp == nshapes) {
                        log.Surv <- slice.shp$log.Surv
                        Ms <- slice.shp$Ms
                        Int <- slice.shp$Int
                        wFun <- slice.shp$wFun
                        wFun2 <- slice.shp$wFun2
                        vl <- slice.shp$vl
                        vls <- slice.shp$vls                    
                    }
                }
            }
        }
        if (control$verbose && !i %% n.batch)
            setTxtProgressBar(pb, i)
        # save results
        if (i %in% resInd) {
            jj <- match(i, resInd)
            if (nbetas1) 
                res.betas1[jj, ] <- betas1
            if (performHC) 
                res.betas2[jj, ] <- betas2            
            if (hasScale)
                res.tau[jj, ] <- tau
            res.b[, , jj] <- b
            if (performHC)
                res.mean_b[, , jj] <- mean_b(indBetas)
            res.invD[jj, ] <- c(invD)
            if (notNullW)
                res.gammas[jj, ] <- gammas
            res.Bs.gammas[jj, ] <- if (rescale_Bs.gammas) {
                tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas
            } else Bs.gammas
            if (baseHazP)
                res.tauBs[jj, ] <- tauBs
            if (estimateAlphas)
                res.alphas[jj, ] <- alphas
            if (paramExtra)
                res.Dalphas[jj, ] <- Dalphas
            if (estimateWeightFun)
                res.shapes[jj, ] <- shapes
            log.pyb <- fastSumID2(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), idFast)
            log.h0s <- if (rescale_Bs.gammas) {
                drop(W2s %*% (tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas))
            } else drop(W2s %*% Bs.gammas)
            if (!paramRE) {
                Mtime <- numeric(n)
                Ms <- numeric(ns)
                if (paramValue) {
                    Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                    Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
                }
                if (paramExtra) {
                    Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                    Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
                }
                if (estimateWeightFun) {
                    wFun <- w * weightFun(u.idGK, shapes, max.time)
                    wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
                    Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
                    Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
                    if (nbetas1) {
                        Xsbetas <- drop(Xs %*% betas1)
                        Xubetas <- drop(Xu %*% betas1)
                        vl <- transFun.value(P * fastSumID2(wFun * (Xsbetas + Zsb), id.GKFast), data.id)
                        vls <- transFun.value(P2 * fastSumID2(wFun2 * (Xubetas + Zub), id.GK2Fast), data.s)
                    } else {
                        vl <- transFun.value(P * fastSumID2(wFun * Zsb, id.GKFast), data.id)
                        vls <- transFun.value(P2 * fastSumID2(wFun2 * Zub, id.GK2Fast), data.s)
                    }
                    Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                    Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
                }                
                if (notNullW) {
                    if (LongFormat) {
                        log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
                    } else {
                        Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                        log.Surv <- expWgammas * Int
                    }
                } else {
                    log.Surv <- Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                }
            } else {
                log.Surv <- Int <- if (notNullW && LongFormat) {
                    P * fastSumID2(w * exp(log.h0s + Wsgammas), id.GKFast)
                } else {
                    P * fastSumID2(w * exp(log.h0s), id.GKFast)
                }
                Mtime <- if (paramSharedRE) {
                    mu_b <- if (performHC) mean_b(indBetas) else betas2
                    drop((b - mu_b) %*% alphas) 
                } else {
                    mu_b <- if (performHC) mean_b(indBetas) else betas2
                    betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
                    if (nbetas1) betas[!indBetasL] <- betas1
                    if (performHC) betas[indBetasL] <- betas2
                    bb_ <- (b - mu_b) + betas[indBetasRE]
                    drop(bb_ %*% alphas)
                }
                log.Surv <- exp(Mtime) * log.Surv
                if (notNullW && !LongFormat)
                    log.Surv <- expWgammas * log.Surv
            }
            log.h <- if (rescale_Bs.gammas) {
                drop(W2 %*% (tchol_CovBs.gammas %*% Bs.gammas + init.Bs.gammas)) + Mtime
            } else drop(W2 %*% Bs.gammas) + Mtime
            if (notNullW)
                log.h <- drop(W %*% gammas) + log.h
            log.ptb <- event * log.h - log.Surv
            mu_b <- if (performHC) mean_b(indBetas) else betas2
            log.pb <- densRE(b, mu = mu_b, invD = invD, log = TRUE, prop = FALSE)
            res.logLik[jj, ] <- log.pyb + log.ptb + log.pb
        }
    })
    if (control$verbose)
        close(pb)
    res.betas <- matrix(0, n.out, if (performHC) nbetas1 + nbetas2 else nbetas1)
    if (nbetas1)
        res.betas[, !indBetasL] <- res.betas1
    if (performHC)
        res.betas[, indBetasL] <- res.betas2
    mcmcOut <- list(betas = res.betas, sigma = if (hasScale) 1/sqrt(res.tau), b = res.b,
                    D = if (ncZ > 1) t(apply(res.invD, 1L, function (x) solve.default(matrix(x, ncZ))))
                    else as.matrix(apply(res.invD, 1L, function (x) solve.default(matrix(x, ncZ)))),
                    gammas = if (notNullW) res.gammas, Bs.gammas = res.Bs.gammas,
                    tauBs = if (baseHazP) res.tauBs,
                    alphas = if (estimateAlphas) res.alphas, Dalphas = if (paramExtra) res.Dalphas,
                    shapes = if (estimateWeightFun) res.shapes)
    mcmcOut <- mcmcOut[!sapply(mcmcOut, is.null)]
    # calculate pD
    D.bar <- - 2 * mean(rowSums(res.logLik, na.rm = TRUE), na.rm = TRUE)
    postMeans <- lapply(mcmcOut, function (x) {
        d <- dim(x)
        if (!is.null(d) && length(d) > 2) apply(x, c(1L, 2L), mean) else colMeans(as.matrix(x))
    })
    dim(postMeans$D) <- c(ncZ, ncZ)
    betas <- postMeans$betas; betas2 <- betas[indBetasL]; betas1 <- betas[!indBetasL]
    if (!performHC) betas2 <- numeric(ncZ)
    sigma <- postMeans$sigma; b <- postMeans$b; D <- postMeans$D
    gammas <- postMeans$gammas; Bs.gammas <- postMeans$Bs.gammas; alphas <- postMeans$alphas
    Dalphas <- postMeans$Dalphas; shapes <- postMeans$shapes
    log.pyb <- fastSumID2(densLong(y.long, drop(X %*% betas1 + rowSums(Z * b[id, , drop = FALSE])),
                                   sigma, log = TRUE, data), idFast)
    log.h0s <- drop(W2s %*% Bs.gammas)
    if (!paramRE) {
        Mtime <- numeric(n)
        Ms <- numeric(ns)
        if (paramValue) {
            if (nbetas1) {
                vl <- transFun.value(drop(Xtime %*% betas1) + rowSums(Ztime * b), data.id)
                vls <- transFun.value(drop(Xs %*% betas1) + rowSums(Zs * b[id.GK, , drop = FALSE]), data.s)
            } else {
                vl <- transFun.value(rowSums(Ztime * b), data.id)
                vls <- transFun.value(rowSums(Zs * b[id.GK, , drop = FALSE]), data.s)                
            }
            Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        if (paramExtra) {
            if (nbetas1) {
                ex <- transFun.extra(drop(Xtime.extra %*% betas1[iF]) +
                                         rowSums(Ztime.extra * b[, iR, drop = FALSE]), data.id)
                exs <- transFun.extra(drop(Xs.extra %*% betas1[iF]) +
                                          rowSums(Zs.extra * b[id.GK, iR, drop = FALSE]), data.s)
            } else {
                ex <- transFun.extra(rowSums(Ztime.extra * b[, iR, drop = FALSE]), data.id)
                exs <- transFun.extra(rowSums(Zs.extra * b[id.GK, iR, drop = FALSE]), data.s)
            }
            Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }
        if (estimateWeightFun) {
            wFun <- w * weightFun(u.idGK, shapes, max.time)
            Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
            wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
            Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
            if (nbetas1) {
                Xsbetas <- drop(Xs %*% betas1)
                Xubetas <- drop(Xu %*% betas1)
                vl <- transFun.value(P * fastSumID2(wFun * (Xsbetas + Zsb), id.GKFast), data.id)
                vls <- transFun.value(P2 * fastSumID2(wFun2 * (Xubetas + Zub), id.GK2Fast), data.s)
            } else {
                vl <- transFun.value(P * fastSumID2(wFun * Zsb, id.GKFast), data.id)
                vls <- transFun.value(P2 * fastSumID2(wFun2 * Zub, id.GK2Fast), data.s)                
            }
            Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        if (notNullW) {
            if (LongFormat) {
                Wsgammas <- drop(Ws %*% gammas)
                log.Surv <- P * fastSumID2(w * exp(log.h0s + Wsgammas + Ms), id.GKFast)
            } else {
                Int <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
                log.Surv <- exp(drop(W %*% gammas)) * Int
            }
        } else {
            log.Surv <- P * fastSumID2(w * exp(log.h0s + Ms), id.GKFast)
        }    
    } else {
        Mtime <- if (paramSharedRE) {
            mu_b <- if (performHC) mean_b(indBetas) else betas2
            drop((b - mu_b) %*% alphas) 
        } else {
            mu_b <- if (performHC) mean_b(indBetas) else betas2
            betas <- numeric(if (performHC) nbetas1 + nbetas2 else nbetas1)
            if (nbetas1) betas[!indBetasL] <- betas1
            if (performHC) betas[indBetasL] <- betas2
            bb_ <- (b - mu_b) + betas[indBetasRE]
            drop(bb_  %*% alphas)
        }
        if (notNullW) {
            if (LongFormat) {
                Wsgammas <- drop(Ws %*% gammas)
                log.Surv <- P * fastSumID2(w * exp(log.h0s + Wsgammas), id.GKFast)
            } else {
                log.Surv <- exp(drop(W %*% gammas)) * P * fastSumID2(w * exp(log.h0s), id.GKFast)
            }
            
        } else {
            log.Surv <- P * fastSumID2(w * exp(log.h0s), id.GKFast)
        }
        log.Surv <- exp(Mtime) * log.Surv
    }
    log.h <- drop(W2 %*% Bs.gammas) + Mtime
    if (notNullW)
        log.h <- drop(W %*% gammas) + log.h
    log.ptb <- event * log.h - log.Surv
    mu_b <- if (performHC) mean_b(indBetas) else rep(0, ncZ)
    log.pb <- densRE(b, mu = mu_b, D = D, log = TRUE, prop = FALSE)
    D.hat <- - 2 * sum(log.pyb + log.ptb + log.pb, na.rm = TRUE)
    pD <- D.bar - D.hat
    indb <- names(mcmcOut) != "b"
    postVarsRE <- apply(res.b, 1L, function (x) var(t(x)))
    dim(postVarsRE)<- c(ncZ, ncZ, n)
    keepD <- length(betas) + 1 + which(!lower.tri(invD, TRUE))
    keepAR <- -seq_len(n.adapt)
    postModes <- lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2L, modes))
    dim(postModes$D) <- c(ncZ, ncZ)
    if (performHC) {
        mcmcOut$b <- res.b - res.mean_b
        postMeans$b <- apply(mcmcOut$b, c(1L, 2L), mean)
    }
    list(mcmc = if (control$keepRE) mcmcOut else mcmcOut[indb], postMeans = postMeans,
         postModes = postModes,
         postVarsRE = postVarsRE,
         StErr = lapply(mcmcOut[indb], stdErr),
         EffectiveSize = lapply(mcmcOut[indb], effectiveSize),
         StDev = lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2L, sd)),
         CIs = lapply(mcmcOut[indb], function (x) 
             apply(as.matrix(x), 2L, quantile, probs = c(0.025, 0.975))),
         Pvalues = lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2L, computeP)),
         vcov = if (ncZ > 1L) var(do.call(cbind, mcmcOut[indb])[, -keepD]) else var(do.call(cbind, mcmcOut[indb])),
         pD = pD, DIC = pD + D.bar, CPO = 1 / colMeans(exp(-res.logLik)),
         LPML = sum(-log(colMeans(exp(-res.logLik))), na.rm = TRUE), time = time,
         scales = list(betas = if (nbetas1) scale.betas1, b = scale.RE, 
                       Bs.gammas = scale.Bs.gammas,
                       gammas = if (notNullW) scale.gammas,
                       alphas = if (estimateAlphas) scale.alphas,
                       Dalphas = if (paramExtra) scale.Dalphas),
         Covs = Covs,
         acceptRates = list(betas = mean(ar.betas[keepAR]), b = colMeans(ar.b[keepAR, ]),
                            D = mean(ar.invD[keepAR]), 
                            Bs.gammas = mean(ar.Bs.gammas[keepAR]),
                            gammas = if (notNullW) mean(ar.gammas[keepAR]),
                            alphas = if (estimateAlphas) mean(ar.alphas[keepAR]),
                            Dalphas = if (paramExtra) mean(ar.Dalphas[keepAR])))
}
