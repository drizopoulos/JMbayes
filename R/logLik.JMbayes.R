logLik.JMbayes <-
function (object, thetas, b, priors = TRUE, marginal.b = TRUE, marginal.thetas = FALSE, 
        full.Laplace = FALSE, useModes = TRUE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("'object' must inherit from class JMbayes.")
    indexRE <- names(object$postMeans) != "b"
    if (missing(thetas)) {
        thetas <- object$postMeans[indexRE]
    }
    if (missing(b))
        b <- ranef(object)
    if (!is.list(thetas) || length(thetas) != length(object$postMeans[indexRE]))
        stop("'thetas' must be a list with the model's parameters with the same structure as ",
            "'object$postMeans'.")    
    if (!is.matrix(b) || dim(b) != dim(ranef(object)))
        stop("'b' must be a numeric matrix with random effects values with the same ",
            "dimensions as 'ranef(object)'.")
    # data and settings
    timeVar <- object$timeVar
    baseHaz <- object$baseHaz
    lag <- object$y$lag
    df.RE <- object$y$df.RE
    param <- object$param
    indFixed <- object$Forms$extraForm$indFixed
    indRandom <- object$Forms$extraForm$indRandom
    id <- object$y$id
    id <- match(id, unique(id))
    id.GK <- object$y$id.GK
    indBetas <- object$y$indBetas
    y <- object$y$y
    Time <- object$y$Time
    event <- object$y$event
    n <- length(Time)
    X <- object$x$X
    Z <- object$x$Z
    W <- object$x$W
    W2 <- object$x$W2
    W2s <- object$x$W2s
    Xtime <- object$x$Xtime
    Ztime <- object$x$Ztime
    Xtime.extra <- object$x$Xtime.extra
    Ztime.extra <- object$x$Ztime.extra
    Xs <- object$x$Xs
    Zs <- object$x$Zs
    Xs.extra <- object$x$Xs.extra
    Zs.extra <- object$x$Zs.extra
    st <- object$x$st
    wk <- object$x$wk
    P <- object$x$P
    wk.long <- rep(wk, n)
    data <- object$Data$data
    data.id <- object$Data$data.id
    data.s <- object$Data$data.s
    densLong <- object$Funs$densLong
    densRE <- object$Funs$densRE
    transFun.value <- object$Funs$transFun.value
    transFun.extra <- object$Funs$transFun.extra
    # parameters
    betas <- thetas$betas
    sigma <- thetas$sigma
    D <- thetas$D
    gammas <- thetas$gammas
    alphas <- thetas$alphas
    Dalphas <- thetas$Dalphas
    Bs.gammas <- thetas$Bs.gammas
    tauBs <- thetas$tauBs
    # log-likelihood
    h <- function (b, individuals = NULL) {
        eta.y <- c(X %*% betas) + rowSums(Z * b[id, , drop = FALSE])
        log.p.y.b <- fastSumID(densLong(y, eta.y, sigma, log = TRUE, data), id)
        eta.t <- if (!is.null(gammas)) c(W %*% gammas) else rep(0, n)
        if (param %in% c("td-value", "td-both")) {
            Y <- transFun.value(c(Xtime %*% betas) + rowSums(Ztime * b), data.id)
            Ys <- transFun.value(c(Xs %*% betas) + rowSums(Zs * b[id.GK, , drop = FALSE]), data.s)
        }
        if (param %in% c("td-extra", "td-both")) {
            Yextra <- transFun.extra(c(Xtime.extra %*% betas[indFixed]) + 
                                         rowSums(Ztime.extra * b[, indRandom, drop = FALSE]), data.id)
            Ys.extra <- transFun.extra(c(Xs.extra %*% betas[indFixed]) + 
                                           rowSums(Zs.extra * b[id.GK, indRandom, drop = FALSE]), data.s)
        }
        longSurv <- c(switch(param,
            "td-value" = as.matrix(Y) %*% alphas, 
            "td-extra" = as.matrix(Yextra) %*% Dalphas,
            "td-both" = as.matrix(Y) %*% alphas + as.matrix(Yextra) %*% Dalphas,
            "shared-betasRE" = (rep(betas[indBetas], each = nrow(b)) + b) %*% alphas,
            "shared-RE" = b %*% alphas))
        longSurv.s <- c(switch(param,
            "td-value" = as.matrix(Ys) %*% alphas, 
            "td-extra" =  as.matrix(Ys.extra) %*% Dalphas,
            "td-both" = as.matrix(Ys) %*% alphas + as.matrix(Ys.extra) %*% Dalphas,
            "shared-betasRE" = c((rep(betas[indBetas], each = nrow(b)) + b) %*% alphas)[id.GK],
            "shared-RE" = c(b %*% alphas)[id.GK]))
        log.hazard <- c(W2 %*% Bs.gammas) + eta.t + longSurv
        log.survival <- - exp(eta.t) * P * fastSumID(wk.long * exp(c(W2s %*% Bs.gammas) + longSurv.s), id.GK)
        log.p.t.b <- event * log.hazard + log.survival
        if (!is.null(individuals))
            log.p.y.b[individuals] + log.p.t.b[individuals]
        else
            log.p.y.b + log.p.t.b
    }
    logLik <- if (!marginal.b) {
        log.p.b <- densRE(b, mu = rep(0, ncol(Z)), D, log = TRUE, prop = FALSE)
        sum(h(b) + log.p.b, na.rm = TRUE)
    } else {
        mean.b <- ranef(object)
        var.b <- attr(ranef(object, postVar = TRUE), "postVar")
        if (full.Laplace) {
            optFun <- function (b, id) {
                b. <- ranef(object)
                b.[id, ] <- b
                log.p.b <- densRE(b, D, log = TRUE, prop = FALSE)
                - h(b., id) - log.p.b[id] 
            }
            for (i in seq_len(n)) {
                opt <- optim(mean.b[i, ], optFun, id = i, method = "BFGS", 
                    hessian = TRUE)
                mean.b[i, ] <- opt$par
                var.b[[i]] <- solve(opt$hessian)
            }
        }
        log.p.b <- densRE(b, mu = rep(0, ncol(Z)), D, log = TRUE, prop = FALSE)
        log.dets.var.b <- apply(var.b, 3, function (x) determinant(x)$modulus)
        sum(0.5 * ncol(b) * log(2 * pi) + 0.5 * log.dets.var.b + 
            h(mean.b) + log.p.b, na.rm = TRUE)
    }
    # priors
    if (priors) {
        priors <- object$priors
        log.betas <- dmvnorm(betas, priors$priorMean.betas, 
            solve(priors$priorTau.betas), log = TRUE)
        log.D <- dwish(solve(D), priors$priorR.invD, priors$priorK.invD, log = TRUE)
        logPrior <- log.betas + log.D
        if (!is.null(sigma)) {
            log.tau <- dgamma(1 / (sigma^2), priors$priorA.tau, 
                              priors$priorB.tau, log = TRUE)
            logPrior <- logPrior + log.tau
        }
        if (!is.null(gammas)) {
            log.gammas <- dmvnorm(gammas, priors$priorMean.gammas, 
                solve(priors$priorTau.gammas), log = TRUE)
            logPrior <- logPrior + log.gammas
        }
        if (!is.null(alphas)) {
            log.alphas <- dmvnorm(alphas, priors$priorMean.alphas, 
                solve(priors$priorTau.alphas), log = TRUE)
            logPrior <- logPrior + log.alphas
        }
        if (!is.null(Dalphas)) {
            log.Dalphas <- dmvnorm(Dalphas, priors$priorMean.Dalphas, 
                solve(priors$priorTau.Dalphas), log = TRUE)
            logPrior <- logPrior + log.Dalphas
        }
        if (baseHaz == "regression-splines") {
            log.Bs.gammas <- dmvnorm(Bs.gammas, priors$priorMean.Bs.gammas, 
                solve(priors$priorTau.Bs.gammas), log = TRUE)
            logPrior <- logPrior + log.Bs.gammas
        } else {
            log.Bs.gammas <- dmvnorm(Bs.gammas, priors$priorMean.Bs.gammas, 
                                     invSigma = tauBs * priors$priorTau.Bs.gammas, 
                                     log = TRUE)
            log.tauBs <- dgamma(tauBs, priors$priorA.tauBs, 
                              priors$priorB.tauBs, log = TRUE)
            logPrior <- logPrior + log.Bs.gammas + log.tauBs
        }
        logLik <- logLik + logPrior
    }
    if (marginal.thetas) {
        tht <- if (useModes) object$postModes else object$postMeans[-3L]
        lL <- logLik(object, thetas = tht)
        var.thetas <- hessian.JMbayes(object, thetas = tht)
        tht$D <- tht$D[lower.tri(tht$D, TRUE)]
        nthetas <- length(unlist(tht))
        logLik <- 0.5 * nthetas * log(2 * pi) - 
            0.5 * determinant(var.thetas)$modulus + lL
    }
    out <- as.vector(logLik)
    attr(out, "df") <- nrow(object$vcov)
    attr(out, "nobs") <- n
    class(out) <- "logLik"
    out
}
