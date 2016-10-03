initSurv_mvJM <- function (Data, control) {
    DF <- data.frame(Time = Data$Time, event = Data$event)
    if (!is.null(W2 <- Data$W2)) {
        DF <- cbind(DF, as.data.frame(W2))
    }
    DF <- cbind(DF, as.data.frame(Data$Wlong))
    fmCox <- coxph(Surv(Time, event) ~ ., data = DF)
    coefs <- coef(fmCox)
    V <- vcov(fmCox)
    out <- NULL
    if (!is.null(W2)) {
        iW <- seq_len(ncol(W2))
        out$gammas <- coefs[iW]
        out$Cov_gammas <- V[iW, iW, drop = FALSE]
        coefs <- coefs[-iW]
        V <- V[-iW, -iW, drop = FALSE]
    }
    iL <- seq_len(ncol(Data$Wlong))
    out$alphas <- coefs[iL]
    out$Cov_alphas <- V[iL, iL, drop = FALSE]
    #############
    W1 <- Data$W1
    W1s <- Data$W1s
    Pw <- Data$Pw
    event <- Data$event
    K <- crossprod(diff(diag(ncol(W1)), diff = control$diff))
    fn1 <- function (Bs_gammas, tauBs) {
        pen <- 0.5 * tauBs * drop(crossprod(Bs_gammas, K %*% Bs_gammas))
        - sum(event * drop(W1 %*% Bs_gammas)) + sum(Pw * exp(W1s %*% Bs_gammas)) + pen
    }
    out$tauBs <- tauBs <- 200
    opt <- optim(rep(0, ncol(W1)), fn1, tauBs = tauBs, method = "BFGS", hessian = TRUE)
    out$Bs_gammas <- opt$par
    out$Cov_Bs_gammas <- solve(opt$hessian)
    out
}
