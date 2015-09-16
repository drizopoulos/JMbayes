marglogLik <-
function (object, newdata, idVar = "id", method = "BFGS", control = NULL) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    timeVar <- object$timeVar
    baseHaz <- object$baseHaz
    df.RE <- object$df.RE
    param <- object$param
    extraForm <- object$Forms$extraForm
    indFixed <- extraForm$indFixed
    indRandom <- extraForm$indRandom
    estimateWeightFun <- object$estimateWeightFun
    weightFun <- object$Funs$weightFun
    max.time <- max(object$y$Time)
    TermsX <- object$Terms$termsYx
    TermsZ <- object$Terms$termsYz
    TermsX.extra <- object$Terms$termsYx.extra
    TermsZ.extra <- object$Terms$termsYz.extra
    mfX <- model.frame(TermsX, data = newdata)
    mfZ <- model.frame(TermsZ, data = newdata)
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$Forms$formYz
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
        rep(TRUE, nrow(newdata))
    } else {
        !seq_len(nrow(newdata)) %in% na.ind
    }
    indBetas <- object$y$indBetas
    id <- as.numeric(unclass(newdata[[idVar]]))
    id <- id. <- match(id, unique(id))
    id <- id[na.ind]
    y <- model.response(mfX)
    X <- model.matrix(formYx, mfX)
    Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
    TermsT <- object$Terms$termsT
    data.id <- newdata[!duplicated(id), ]
    data.s <- data.id[rep(1:nrow(data.id), each = 15), ]
    idT <- data.id[[idVar]]
    idT <- match(idT, unique(idT))
    ids <- data.s[[idVar]]
    ids <- match(ids, unique(ids))
    mfT <- model.frame(delete.response(TermsT), data = data.id)
    formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
        strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else {
        tt <- attr(delete.response(TermsT), "term.labels")
        if (length(tt)) reformulate(tt) else reformulate("1")
    }
    W <- model.matrix(formT, mfT)[, -1, drop = FALSE]
    last.time <- tapply(newdata[[timeVar]], id., tail, n = 1)
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    lag <- object$y$lag
    densLong <- object$Funs$densLong
    densRE <- object$Funs$densRE
    transFun.value <- object$Funs$transFun.value
    transFun.extra <- object$Funs$transFun.extra
    # list of parameters
    list.thetas <- c(object$postModes, list(ranef = rep(0, ncz) ))
    if (!is.null(list.thetas$sigma))
        list.thetas$sigma <- log(list.thetas$sigma)
    if (!is.null(list.thetas$tauBs))
        list.thetas$tauBs <- log(list.thetas$tauBs)
    list.thetas$D <- chol.transf(list.thetas$D)
    thetas.b <- unlist(as.relistable(list.thetas))
    # construct model matrices to calculate the log.posterior
    environment(ModelMats) <- environment()
    survMats.last <- vector("list", n.tp)
    for (i in seq_along(last.time)) {
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i)
    }
    logPost <- function (thetas.b, ii, transform = TRUE) {
        # parameters
        tht <- relist(thetas.b, list.thetas)
        betas <- tht$betas        
        gammas <- tht$gammas
        alphas <- tht$alphas
        Dalphas <- tht$Dalphas
        Bs.gammas <- tht$Bs.gammas
        sigma <- tht$sigma
        tauBs <- tht$tauBs
        b <- tht$ranef
        if (transform) {
            if (!is.null(sigma)) sigma <- exp(tht$sigma)
            if (!is.null(tauBs)) tauBs <- exp(tauBs)
            D <- chol.transf(tht$D)
        } else {
            sigma <- tht$sigma
            D <- matrix(0, ncz, ncz)
            D[lower.tri(D, TRUE)] <- tht$D
            D <- D + t(D)
            diag(D) <- diag(D) / 2
        }
        # log-likelihood contributions
        id.i <- id %in% ii
        idT.i <- idT %in% ii
        ids.i <- ids %in% ii
        X.i <- X[id.i, , drop = FALSE]
        Z.i <- Z[id.i, , drop = FALSE]
        mu.y <- as.vector(X.i %*% betas + Z.i %*% b)
        logY <- densLong(y[id.i], mu.y, sigma, log = TRUE)
        log.p.yb <- sum(logY)
        log.p.b <- densRE(b, mu = rep(0, ncz), D = D, log = TRUE, prop = FALSE)
        st <- survMats.last[[ii]]$st
        wk <- survMats.last[[ii]]$wk
        P <- survMats.last[[ii]]$P
        Xs <- survMats.last[[ii]]$Xs
        Zs <- survMats.last[[ii]]$Zs
        Xs.extra <- survMats.last[[ii]]$Xs.extra
        Zs.extra <- survMats.last[[ii]]$Zs.extra
        W2s <- survMats.last[[ii]]$W2s
        ind <- survMats.last[[ii]]$ind
        if (param %in% c("td-value", "td-both"))
            Ys <- transFun.value(as.vector(Xs %*% betas + Zs %*% b), data.s[ids.i, ])
        if (param %in% c("td-extra", "td-both"))
            Ys.extra <- transFun.extra(as.vector(Xs.extra %*% betas[indFixed] + 
                                         Zs.extra %*% b[indRandom]), data.s[ids.i, ])
        tt <- switch(param,
            "td-value" = as.matrix(Ys) %*% alphas, 
            "td-extra" = as.matrix(Ys.extra) %*% Dalphas,
            "td-both" = as.matrix(Ys) %*% alphas + as.matrix(Ys.extra) %*% Dalphas,
            "shared-betasRE" = rep(sum((betas[indBetas] + b) * alphas), length(st)),        
            "shared-RE" = rep(sum(b * alphas), length(st)))
        eta.tw <- if (ncol(W)) {
                as.vector(W[ii, , drop = FALSE] %*% gammas)
        } else 0
        Vi <- exp(c(W2s %*% Bs.gammas) + tt)
        idT <- rep(seq_along(P), each = 15)
        log.survival <- - sum(exp(eta.tw) * P * tapply(wk * Vi, idT, sum))
        if (all(st == 0))
            log.survival <- 1
        logLik <- log.p.yb + log.survival + log.p.b
        # Priors
        priors <- object$priors
        log.betas <- dmvnorm(betas, priors$priorMean.betas, 
            solve(priors$priorTau.betas), log = TRUE)
        log.tau <- dgamma(1 / (sigma^2), priors$priorA.tau, 
            priors$priorB.tau, log = TRUE)
        log.D <- dwish(solve(D), priors$priorR.invD, priors$priorK.invD, log = TRUE)
        logPrior <- log.betas + log.tau + log.D
        if (!is.null(gammas)) {
            ind <- colSums(object$x$W == 0) == nrow(object$x$W)
            log.gammas <- dmvnorm(gammas, priors$priorMean.gammas[!ind], 
                solve(priors$priorTau.gammas)[!ind, !ind], log = TRUE)
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
        - as.vector(logLik + logPrior)
    }
    score.logPost <- function (thetas.b, ii, transform = TRUE) {
        fd(thetas.b, logPost, ii = ii, transform = transform)
    }
    w <- numeric(length(last.time))
    con <- list(maxit = 200, parscale = rep(0.001, length(thetas.b)))
    con[names(control)] <- control
    for (i in seq_along(last.time)) {
        w[i] <- tryCatch({
            opt <- optim(thetas.b, logPost, score.logPost, ii = i, transform = TRUE, 
                method = method, control = con, hessian = TRUE)
            opt.thetas <- relist(opt$par, list.thetas)
            opt.thetas$sigma <- exp(opt.thetas$sigma)
            opt.thetas$D <- chol.transf(opt.thetas$D)
            opt.thetas$D <- opt.thetas$D[lower.tri(opt.thetas$D, TRUE)]
            opt.thetas <- unlist(as.relistable(opt.thetas))
            H <- opt$hessian
            as.vector(0.5 * length(unlist(list.thetas)) * log(2 * pi) - 
                0.5 * determinant(H)$modulus - opt$value)
        }, error = function (e) NA)
    }
    w
}
