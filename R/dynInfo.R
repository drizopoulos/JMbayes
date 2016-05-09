dynInfo <- function (object, newdata, Dt, K = 5, M = 500, idVar = "id", 
                     simulateFun = function (eta, scale) rnorm(length(eta), eta, scale), 
                     seed = 1L) {
    if (!inherits(object, "JMbayes"))
        stop("'object' must inherit from class JMbayes.")
    set.seed(seed)
    # extract components from fitted joint model
    timeVar <- object$timeVar
    df.RE <- object$y$df.RE
    param <- object$param
    densLong <- object$Funs$densLong
    hasScale <- object$Funs$hasScale
    densRE <- object$Funs$densRE
    transFun.value <- object$Funs$transFun.value
    transFun.extra <- object$Funs$transFun.extra
    extraForm <- object$Forms$extraForm
    indFixed <- extraForm$indFixed
    indRandom <- extraForm$indRandom
    performHC <- object$control$performHC
    TermsX <- object$Terms$termsYx
    TermsZ <- object$Terms$termsYz
    TermsX.extra <- object$Terms$termsYx.extra
    TermsZ.extra <- object$Terms$termsYz.extra
    mfX <- model.frame.default(TermsX, data = newdata)
    mfZ <- model.frame.default(TermsZ, data = newdata)
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$Forms$formYz
    estimateWeightFun <- object$estimateWeightFun
    weightFun <- object$Funs$weightFun
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
        rep(TRUE, nrow(newdata))
    } else {
        !seq_len(nrow(newdata)) %in% na.ind
    }
    id <- as.numeric(unclass(newdata[[idVar]]))
    id <- match(id, unique(id))
    id <- id. <- id[na.ind]
    y <- yy <- model.response(mfX)
    X <- XX <-  model.matrix.default(formYx, mfX)
    Z <- ZZ <- model.matrix.default(formYz, mfZ)[na.ind, , drop = FALSE]
    TermsT <- object$Terms$termsT
    # extract time and variables from newdata
    TimeVar <- all.vars(TermsT)[1L]
    eventVar <- all.vars(TermsT)[2L]
    respVar <- all.vars(TermsX)[1L]
    maxTime <- max(object$y$Time)
    max_time <- max(newdata[[timeVar]])
    times <- seq(max_time, max_time + Dt, len = K + 1)[-1L]
    ntimes <- length(times)
    data.id <- newdata[tapply(row.names(newdata), id, tail, n = 1L),]
    data.s <- data.id[rep(1:nrow(data.id), each = object$control$GQsurv.k), ]
    idT <- data.id[[idVar]]
    idT <- match(idT, unique(idT))
    ids <- data.s[[idVar]]
    ids <- match(ids, unique(ids))
    data.p <- data.id[rep(1:nrow(data.id), each = ntimes), ]
    data.p[[timeVar]] <- times
    mfXpred <- model.frame.default(TermsX, data = data.p)
    mfZpred <- model.frame.default(TermsZ, data = data.p)
    Xpred <- model.matrix.default(formYx, mfXpred)
    Zpred <- model.matrix.default(formYz, mfZpred)
    mfT <- model.frame.default(delete.response(TermsT), data = data.id)
    tt <- attr(delete.response(TermsT), "term.labels")
    formT <- if (length(tt)) reformulate(tt) else reformulate("1")
    W <- model.matrix.default(formT, mfT)[, -1L, drop = FALSE]
    ########
    n <- nrow(data.id)
    n.tp <- length(times)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    if (ncww == 0L)
        W <- NULL
    lag <- object$y$lag
    betas <- object$postMeans$betas
    sigma <- object$postMeans$sigma
    D <- object$postMeans$D
    gammas <- object$postMeans$gammas
    alphas <- object$postMeans$alphas
    Dalphas <- object$postMeans$Dalphas
    shapes <- object$postMeans$shapes
    Bs.gammas <- object$postMeans$Bs.gammas
    list.thetas <- list(betas = betas, sigma = sigma, gammas = gammas, alphas = alphas, 
                        Dalphas = Dalphas, shapes = shapes, Bs.gammas = Bs.gammas, D = D)
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    ########
    environment(log.posterior.b) <- environment(S.b) <- environment(logh.b) <- environment()
    environment(hMats) <- environment(ModelMats) <- environment()
    obs.times <- split(newdata[[timeVar]][na.ind], id)
    survMats <- lapply(c(max_time, times), ModelMats, ii = 1)
    u_times <- lapply(lapply(times, seq, to = maxTime * 1.1, length.out = 31), tail, n = -1)
    #hh <- function (t) lapply(t, ModelMats, ii = 1)
    #survMats_samp <- lapply(u_times, hh)
    #hazMats <- lapply(times, hMats)
    # calculate the Empirical Bayes estimates and their (scaled) variance
    betas.new <- betas
    sigma.new <- sigma
    D.new <- D
    gammas.new <- gammas
    alphas.new <- alphas
    Dalphas.new <- Dalphas
    shapes.new <- shapes
    Bs.gammas.new <- Bs.gammas
    ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y, Mats = tt, ii = i)
    start <- rep(0, ncz)
    opt <- try(optim(start, ff, y = y, tt = survMats, i = 1, 
                     method = "BFGS", hessian = TRUE), silent = TRUE)
    if (inherits(opt, "try-error")) {
        gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, i = i)
        opt <- optim(start, ff, gg, y = y, tt = survMats, 
                     i = 1, method = "BFGS", hessian = TRUE, 
                     control = list(parscale = rep(0.1, ncz)))
    } 
    modes.b <- opt$par
    invVars.b <- opt$hessian / 1.5
    Vars.b <- 1.5 * solve(opt$hessian)
    b.old <- b.new <- modes.b
    b.old1 <- b.new1 <- modes.b
    b.old2 <- b.new2 <- modes.b
    mcmc <- object$mcmc
    mcmc <- mcmc[names(mcmc) != "b"]
    samples <- sample(nrow(mcmc$betas), 3 * M, replace = TRUE)
    mcmc[] <- lapply(mcmc, function (x) x[samples, , drop = FALSE])
    proposed.b <- rmvt(n = M, df = 4, mu = modes.b, Sigma = Vars.b)
    dmvt.proposed <- dmvt(proposed.b, mu = modes.b, Sigma = Vars.b, df = 4, log = TRUE)
    proposed.b1 <- rmvt(n = M, df = 4, mu = modes.b, Sigma = Vars.b)
    dmvt.proposed1 <- dmvt(proposed.b1, mu = modes.b, Sigma = Vars.b, df = 4, log = TRUE)
    proposed.b2 <- rmvt(n = M, df = 4, mu = modes.b, Sigma = Vars.b)
    dmvt.proposed2 <- dmvt(proposed.b2, mu = modes.b, Sigma = Vars.b, df = 4, log = TRUE)    
    # loop over time points
    log.p_Tj <- function (Tj) {
        log.S_ti <- S.b(times[1], b.new2, i = 1, survMats[[1]], log = TRUE)
        log.S_Tj <- S.b(Tj, b.new2, i = 1, ModelMats(Tj, 1), log = TRUE)
        log.h_Tj <- logh.b(b.new2, hMats(Tj))
        log.h_Tj + log.S_Tj - log.S_ti
    }
    sfit <- survfitJM(object, newdata = newdata, M = M, init.b = rbind(modes.b),
                      survTimes = times, idVar = idVar)
    sfit <- 1 - as.vector(sfit$summaries[[1]][, "Mean"])
    info.times <- matrix(0, M, ntimes)
    for (ti in seq_len(ntimes)) {
        # Monte Carlo scheme
        old_Tj <- Tj <- 1.1 * times[ti]
        count <- count.b <- count.b1 <- count.b2 <- 0.0
        info <- numeric(M)
        for (m in seq_len(M)) {
            # Step 1-1: Simulate parameter values from [theta | D_n]
            betas.new <- mcmc$betas[m, ]
            if (hasScale)
                sigma.new <- mcmc$sigma[m, ]
            if (!is.null(W))
                gammas.new <- mcmc$gammas[m, ]
            if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
                alphas.new <- mcmc$alpha[m, ]
            if (param %in% c("td-extra", "td-both"))
                Dalphas.new <- mcmc$Dalphas[m, ]
            if (estimateWeightFun)
                shapes.new <- mcmc$shapes[m, ]
            D.new <- mcmc$D[m, ]; dim(D.new) <- dim(D)
            Bs.gammas.new <- mcmc$Bs.gammas[m, ]
            # Step 1-2: Simulate b from [b | T_j > t, Y_j(t)]
            id <- id.
            y <- yy
            X <- XX
            Z <- ZZ
            p.b <- proposed.b[m, ]
            dmvt.old <- dmvt(b.old, modes.b, invSigma = invVars.b, df = 4, log = TRUE)
            dmvt.prop <- dmvt.proposed[m]
            a <- min(exp(log.posterior.b(p.b, y, list(survMats[[1]]), ii = 1) + dmvt.old - 
                             log.posterior.b(b.old, y, list(survMats[[1]]), ii = 1) - dmvt.prop), 1)
            ind <- runif(1) <= a
            if (!is.na(ind) && ind) {
                b.new <- p.b
                count.b <- count.b + 1
            }
            b.old <- b.new
            # Step 1-3: Simulate y_j(u)
            eta.y <- drop(Xpred %*% betas.new + Zpred %*% b.new)[ti]
            y.new <- simulateFun(eta.y, sigma.new)
            ##################################################################################
            # Step 2-1: Simulate parameter values from [theta | D_n]
            betas.new <- mcmc$betas[m + M, ]
            if (hasScale)
                sigma.new <- mcmc$sigma[m + M, ]
            if (!is.null(W))
                gammas.new <- mcmc$gammas[m + M, ]
            if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
                alphas.new <- mcmc$alpha[m + M, ]
            if (param %in% c("td-extra", "td-both"))
                Dalphas.new <- mcmc$Dalphas[m + M, ]
            if (estimateWeightFun)
                shapes.new <- mcmc$shapes[m + M, ]
            D.new <- mcmc$D[m + M, ]; dim(D.new) <- dim(D)
            Bs.gammas.new <- mcmc$Bs.gammas[m + M, ]
            # Step 2-2: Simulate b from [b | T_j > t, {Y_j(t), y_j(u)}]
            id <- c(id., tail(id, 1))
            y <- c(yy, y.new)
            X <- rbind(XX, Xpred[ti, ])
            Z <- rbind(ZZ, Zpred[ti, ])
            p.b1 <- proposed.b1[m, ]
            dmvt.old1 <- dmvt(b.old1, modes.b, invSigma = invVars.b, df = 4, log = TRUE)
            dmvt.prop1 <- dmvt.proposed1[m]
            a1 <- min(exp(log.posterior.b(p.b1, y, list(survMats[[1]]), ii = 1) + dmvt.old1 - 
                              log.posterior.b(b.old1, y, list(survMats[[1]]), ii = 1) - dmvt.prop1), 1)
            ind1 <- runif(1) <= a1
            if (!is.na(ind1) && ind1) {
                b.new1 <- p.b1
                count.b1 <- count.b1 + 1
            }
            b.old1 <- b.new1
            # Step 2-3: Simulate T_j^* from [T_j^* | T_j > u, {Y_j(t), y_j(u)}]
            prop_Tj <- runif(1, times[ti], maxTime * 1.1)
            aa <- min(exp(log.p_Tj(prop_Tj) - log.p_Tj(old_Tj)), 1)
            ind <- runif(1) <= aa
            if (!is.na(ind) && ind) {
                Tj <- prop_Tj
                count <- count + 1
            }
            old_Tj <- Tj
            ##################################################################################
            # Step 3-1: Simulate parameter values from [theta | D_n]
            betas.new <- mcmc$betas[m + 2*M, ]
            if (hasScale)
                sigma.new <- mcmc$sigma[m + 2*M, ]
            if (!is.null(W))
                gammas.new <- mcmc$gammas[m + 2*M, ]
            if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
                alphas.new <- mcmc$alpha[m + 2*M, ]
            if (param %in% c("td-extra", "td-both"))
                Dalphas.new <- mcmc$Dalphas[m + 2*M, ]
            if (estimateWeightFun)
                shapes.new <- mcmc$shapes[m + 2*M, ]
            D.new <- mcmc$D[m + 2*M, ]; dim(D.new) <- dim(D)
            Bs.gammas.new <- mcmc$Bs.gammas[m +2*M, ]
            # Step 3-2: Simulate b from [b | T_j > t, {Y_j(t), y_j(u)}]
            p.b2 <- proposed.b2[m, ]
            dmvt.old2 <- dmvt(b.old2, modes.b, invSigma = invVars.b, df = 4, log = TRUE)
            dmvt.prop2 <- dmvt.proposed2[m]
            a2 <- min(exp(log.posterior.b(p.b2, y, list(survMats[[1 + ti]]), ii = 1) + dmvt.old2 - 
                              log.posterior.b(b.old2, y, list(survMats[[1 + ti]]), ii = 1) - dmvt.prop2), 1)
            ind2 <- runif(1) <= a2
            if (!is.na(ind2) && ind2) {
                b.new2 <- p.b2
                count.b2 <- count.b2 + 1
            }
            b.old2 <- b.new2
            # Step 3-3: Calculate p(T_j | T_j > t, b)
            log.S_ti <- S.b(times[ti], b.new2, i = 1, survMats[[1 + ti]], log = TRUE)
            log.S_Tj <- S.b(Tj, b.new2, i = 1, ModelMats(Tj, 1), log = TRUE)
            log.h_Tj <- logh.b(b.new2, hMats(Tj))
            info[m] <- log.h_Tj + log.S_Tj - log.S_ti
        }
        #print(count/M)
        #print(count.b/M)
        #print(count.b1/M)
        #print(count.b2/M)
        info.times[, ti] <- info
    }
    infoSum <- apply(info.times, 2, median, na.rm = TRUE)
    stars <- rep(" ", length.out = ntimes)
    stars[which.max(infoSum)] <- "*"
    d <- data.frame(times = times, Info = infoSum, pi = sfit, 
                    " " = stars, check.names = FALSE)
    rm(list = ".Random.seed", envir = globalenv())
    list(summary = d, full.results = info.times)
}
