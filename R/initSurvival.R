initSurvival <- function (Time, event, id, W2, W2s, P, wk, id.GK, times, b = NULL, betas = NULL, 
                          indBetas = NULL, W = NULL, baseHaz = NULL, diff = NULL,
                          Data = NULL, param = NULL, long = NULL, long.extra = NULL, 
                          transFun.value = NULL, transFun.extra = NULL,
                          vl = NULL, vls = NULL, ex = NULL, exs = NULL) {
    nTime <- length(Time)
    nLong <- length(id)
    if (param == "shared-betasRE" || param == "shared-RE") {
        DF <- data.frame(Time = Time, event = event)
        Wdat <- as.data.frame(W)
        if (!is.null(W))
            DF <- cbind(DF, Wdat)
        DF <- cbind(DF, b)
        tdCox <- coxph(Surv(Time, event) ~ ., data = DF)        
    } else {
        DF <- data.frame(id = id, Time = Time[id], event = event[id])
        if (!is.null(W)) {
            Wdat <- as.data.frame(W)
            DF <- cbind(DF, Wdat[id, ])
        }
        if (!is.null(long)) {
            long <- as.data.frame(transFun.value(long, Data$data))
            DF <- cbind(DF, long)
        }
        if (!is.null(long.extra)) {
            long.extra <- as.data.frame(transFun.extra(long.extra, Data$data))
            DF <- cbind(DF, long.extra)
        }
        row.names(DF) <- 1:nLong
        DF$start <- times
        splitID <- split(DF[c("start", "Time")], DF$id)
        DF$stop <- unlist(lapply(splitID, function (d) c(d$start[-1], d$Time[1])))
        DF$event <- with(DF, ave(event, id, FUN = function (x) c(rep(0, length(x)-1), x[1])))
        DF <- DF[!names(DF) %in% c("Time", "id")]
        tdCox <- coxph(Surv(start, stop, event) ~ ., data = DF[DF$stop > DF$start, ])
    }
    coefs <- coef(tdCox)
    V <- vcov(tdCox)
    out <- NULL
    if (!is.null(W)) {
        iW <- 1:ncol(W)
        out$gammas <- coefs[iW]
        out$cov.gammas <- V[iW, iW]
        coefs <- coefs[-iW]
        V <- V[-iW, -iW, drop = FALSE]
    }
    if (!is.null(long) || param %in% c("shared-betasRE", "shared-RE")) {
        iL <- if (!is.null(long)) 1:ncol(long) else 1:ncol(b)
        out$alphas <- coefs[iL]
        out$cov.alphas <- V[iL, iL]
        coefs <- coefs[-iL]
        V <- V[-iL, -iL, drop = FALSE]
    }
    if (!is.null(long.extra)) {
        iLe <- 1:ncol(long.extra)
        out$Dalphas <- coefs[iLe]
        out$cov.Dalphas <- V[iLe, iLe]
    }
    if (baseHaz == "P-splines") {
        K <- crossprod(diff(diag(ncol(W2)), diff = diff))
        fn1 <- function (Bs.gammas, tauBs) {
            pen <- 0.5 * tauBs * drop(crossprod(Bs.gammas, K %*% Bs.gammas))
            - sum(event * drop(W2 %*% Bs.gammas) - P * fastSumID(rep(wk, nTime) * exp(drop(W2s %*% Bs.gammas)), id.GK)) + pen
        }
        out$tauBs <- tauBs <- 200
        opt <- optim(rep(0, ncol(W2)), fn1, tauBs = tauBs, method = "BFGS", hessian = TRUE)
    } else {
        fn2 <- function (Bs.gammas) {
            - sum(event * drop(W2 %*% Bs.gammas) - P * fastSumID(rep(wk, nTime) * exp(drop(W2s %*% Bs.gammas)), id.GK))
        }
        opt <- optim(rep(0, ncol(W2)), fn2, method = "BFGS", hessian = TRUE)        
    }
    out$Bs.gammas <- opt$par
    out$cov.Bs.gammas <- solve(opt$hessian)
    ind <- !names(out) %in% c("cov.gammas", "cov.alphas", "cov.Dalphas", "cov.Bs.gammas", "tauBs")
    out.vec <- unlist(as.relistable(out[ind]))
    fn3 <- function (thetas) {
        thetas <- relist(thetas, out[ind])
        gammas <- thetas$gammas
        Bs.gammas <- thetas$Bs.gammas
        alphas <- thetas$alphas
        Dalphas <- thetas$Dalphas
        ####
        ns <- length(id.GK)
        Mtime <- rep(0, nTime)
        Ms <- rep(0, ns)
        if (param %in% c("td-value", "td-both")) {
            Mtime <- Mtime + if (is.matrix(vl)) c(vl %*% alphas) else vl * alphas
            Ms <- Ms + if (is.matrix(vls)) c(vls %*% alphas) else vls * alphas
        }
        if (param %in% c("td-extra", "td-both")) {
            Mtime <- Mtime + if (is.matrix(ex)) c(ex %*% Dalphas) else ex * Dalphas
            Ms <- Ms + if (is.matrix(exs)) c(exs %*% Dalphas) else exs * Dalphas
        }
        if (param == "shared-RE")
            Mtime <- c(b %*% alphas)
        if (param == "shared-betasRE")
            Mtime <- c((rep(betas[indBetas], each = nrow(b)) + b) %*% alphas)        
        log.h0 <- c(W2 %*% Bs.gammas)
        log.h0s <- c(W2s %*% Bs.gammas)
        etaW <- if (is.null(W)) rep(0, nTime) else c(W %*% gammas)
        log.h <- log.h0 + etaW + Mtime
        if (param == "shared-betasRE" || param == "shared-RE")
            etaW <- etaW + Mtime
        log.S <- exp(etaW) * P * fastSumID(rep(wk, nTime) * exp(log.h0s + Ms), id.GK)
        pen <- if (baseHaz == "P-splines") {
            0.5 * tauBs * c(crossprod(Bs.gammas, K %*% Bs.gammas))
        } else 0
        - sum(event * log.h - log.S, na.rm = TRUE) + pen
    }
    test <- try(opt2 <- suppressWarnings(optim(out.vec, fn3, method = "BFGS", hessian = TRUE, 
                                   control = list(parscale = rep(0.1, length(out.vec))))), silent = TRUE)
    if (!inherits(test, "try-error") && !opt2$convergence && eigen(opt2$hessian, TRUE)$values > 0) {
        res <- relist(opt2$par, out[ind])
        out[names(res)] <- res
        V <- solve(nearPD(opt2$hessian))
        if (!is.null(W)) {
            iW <- 1:ncol(W)
            out$cov.gammas <- V[iW, iW]
            V <- V[-iW, -iW, drop = FALSE]
        }
        if (param %in% c("td-value", "td-both")) {
            iL <- 1:ncol(long)
            out$cov.alphas <- V[iL, iL]
            V <- V[-iL, -iL, drop = FALSE]
        }
        if (param %in% c("td-extra", "td-both")) {
            iLe <- 1:ncol(long.extra)
            out$cov.Dalphas <- V[iLe, iLe]
            V <- V[-iLe, -iLe, drop = FALSE]
        }
        if (param %in% c("shared-betasRE", "shared-RE")) {
            iSRE <- 1:ncol(b)
            out$cov.alphas <- V[iSRE, iSRE]
            V <- V[-iSRE, -iSRE, drop = FALSE]
        }
        out$cov.Bs.gammas <- V
    }
    if (!is.null(out$alphas) && is.na(out$alphas)) {
        out$alphas <- rep(0, length(out$alphas))
        out$cov.alphas <- diag(0.1, length(out$alphas))
    }
    if (!is.null(out$Dalphas) && is.na(out$Dalphas)) {
        out$Dalphas <- rep(0, length(out$Dalphas))
        out$cov.Dalphas <- diag(0.1, length(out$Dalphas))
    }
    out
}
