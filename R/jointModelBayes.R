jointModelBayes <- function (lmeObject, survObject, timeVar,
                             param = c("td-value", "td-extra", "td-both", "shared-betasRE", "shared-RE"),
                             extraForm = NULL, baseHaz = c("P-splines", "regression-splines"), 
                             transFun = NULL, densLong = NULL, lag = 0, df.RE = NULL,
                             estimateWeightFun = FALSE, weightFun = NULL, init = NULL, priors = NULL, scales = NULL, 
                             control = list(), ...) {
    cl <- match.call()
    param <- match.arg(param)
    baseHaz <- match.arg(baseHaz)
    if (!inherits(lmeObject, "lme"))
        stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1)
        stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct))
        warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!inherits(survObject, "coxph"))
        stop("\n'survObject' must inherit from class coxph.")
    if (is.null(survObject$x))
        stop("\nuse argument 'x = TRUE' in coxph().")
    if (length(timeVar) != 1L || !is.character(timeVar))
        stop("\n'timeVar' must be a character string.")
    if (param %in% c("td-extra", "td-both") && is.null(extraForm)) {
        stop("\nwhen parameterization is 'td-extra' or 'td-both' you need to specify the 'extraForm' argument.")
    }
    if (!param %in% c("td-extra", "td-both") && !is.null(extraForm)) {
        stop("\nyou have defined 'extraForm' but the parameterization is neither 'td-extra' nor 'td-both'.")
    }
    if (param %in% c("td-extra", "td-both") && !is.list(extraForm)) {
        stop("\nthe 'extraForm' argument must be a list with components 'fixed' (a formula),\n\t'indFixed'",
             "(a numeric vector), 'random' (a formula) and 'indRandom' (a numeric vector).")
    }
    # extract response & design matrix survival process
    formT <- formula(survObject)
    W <- survObject$x
    if (!length(W))
        W <- NULL
    SurvInf <- survObject$y
    typeSurvInf <- attr(SurvInf, "type")
    if (typeSurvInf == "right") {
        Time <- SurvInf[, "time"]
        Time[Time < 1e-04] <- 1e-04
        nT <- length(Time)
        event <- SurvInf[, "status"]
        LongFormat <- FALSE
    }
    if (typeSurvInf == "counting") {
        if (is.null(survObject$model))
            stop("\nplease refit the Cox model including in the ", 
                 "call to coxph() the argument 'model = TRUE'.")
        idT <- if (!is.null(survObject$model$cluster)) {
            as.vector(unclass(survObject$model$cluster))
        } else {
            seq_len(nrow(survObject$model))
        }
        strata <- if (!is.null(survObject$model$strata)) {
            as.vector(unclass(survObject$model$strata))
        } else {
            seq_len(nrow(survObject$model))
        }
        idT <- match(idT, unique(idT))
        LongFormat <- length(idT) > length(unique(idT))
        TimeL <- TimeLl <- SurvInf[, "start"]
        TimeL <- tapply(TimeL, idT, head, n = 1)
        anyLeftTrunc <- any(TimeL > 1e-07)
        TimeR <- SurvInf[, "stop"]
        TimeR[TimeR < 1e-04] <- 1e-04
        Time <- tapply(TimeR, idT, tail, n = 1)
        nT <- length(Time)
        eventLong <- SurvInf[, "status"]
        event  <- tapply(eventLong, idT, tail, n = 1)
    }
    # longitudinal process
    idOrig <- lmeObject$groups[[1L]]
    id <- as.vector(unclass(idOrig))
    b <- data.matrix(ranef(lmeObject))
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT)
        stop("sample sizes in the longitudinal and event processes differ.\n")
    data <- lmeObject$data
    if (!timeVar %in% names(data))
        stop("\n'timeVar' does not correspond to one of the columns in the model.frame of the 'lmeObject'.")
    times <- data[[timeVar]]
    # check if there are any longitudinal measurements after the event times
    max.timeY <- tapply(times, id, max)
    if (!all(Time >= max.timeY)) {
        idnams <- factor(idOrig)
        stop("\nit seems that there are longitudinal measurements taken after the event times for some subjects ",
             "(i.e., check subject(s): ", paste(levels(idnams)[(Time < max.timeY)], collapse = ", "), ").")
    }
    formYx <- formula(lmeObject)
    mfX <- model.frame(terms(formYx), data = data)
    TermsX <- attr(mfX, "terms")
    X <- model.matrix(formYx, mfX)
    offset <- model.offset(mfX)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ)
    y.long <- model.response(mfX, "numeric")
    long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ])    
    if (densLongCheck <- is.null(densLong)) {
        if (is.null(lmeObject$family)) {
            densLong <- function (y, eta.y, scale, log = FALSE, data) {
                dnorm(x = y, mean = eta.y, sd = scale, log = log)
            }
        } else {
            stop("you should define 'densLong' appropriately.\n")
        }
    } else {
        if (!is.function(densLong) || 
                !names(formals(densLong)) %in% c("y", "eta.y", "scale", "log", "data"))
            stop("invalid specification of densLong(); check the help file.")
    }
    # density of the random effects
    densRE <- if (is.null(df.RE)) {
        function (b, mu, D = NULL, invD = NULL, log = FALSE, prop = TRUE) {
            dmvnorm(b, mu = mu, Sigma = D, invSigma = invD, log = log,
                    prop = prop)
        }
    } else {
        ddRE <- function (b, mu = NULL, D = NULL, invD = NULL, log = FALSE, df, prop = TRUE) {
            if (is.null(mu)) {
                nn <- if (is.matrix(b)) ncol(b) else length(b)
                mu <- numeric(nn)
            }
            dmvt(b, mu = mu, Sigma = D, invSigma = invD, df = df, log = log, 
                 prop = prop)
        }
        formals(ddRE)$"df" <- df.RE
        ddRE
    }
    # posterior variances random effects
    D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
                lmeObject$sigma^2)[[1]]
    invD <- solve(D)
    sigma <- lmeObject$sigma
    sigma2 <- sigma * sigma    
    if (densLongCheck && is.null(df.RE)) {
        Cov.postRE <- vector("list", nY)
        for (i in seq_len(nY)) {
            Z.i <- Z[id == i, , drop = FALSE]
            Cov.postRE[[i]] <- solve.default(crossprod(Z.i) / sigma2 + invD)
        }
    } else {
        betas <- fixef(lmeObject)
        logpostRE <- function (b) {
            eta.yi <- drop(X.i %*% betas) + drop(Z.i %*% b)
            dY <- densLong(y.i, eta.yi, scale = sigma, log = TRUE, data[id.i, , drop = FALSE])
            dRE <- densRE(b, mu = rep(0, ncol(Z.i)), invD = invD, log = TRUE)
            - sum(dY, dRE, na.rm = TRUE)
        }
        Cov.postRE <- vector("list", nY)
        for (i in seq_len(nY)) {
            id.i <- id == i
            y.i <- y.long[id.i]
            X.i <- X[id.i, , drop = FALSE]
            Z.i <- Z[id.i, , drop = FALSE]
            opt <- optim(b[i, ], logpostRE, method = "BFGS", hessian = TRUE)
            b[i, ] <- opt$par
            Cov.postRE[[i]] <- solve(opt$hessian)
        }
    }
    # check parameterization and hierarchical centering
    check_names <- all(colnames(Z) %in% colnames(X))
    if (check_names) {
        performHC <- max(diag(D)) > sigma2 / nY && is.null(df.RE)
        has_interceptX <- attr(TermsX, "intercept")
        has_interceptZ <- attr(TermsZ, "intercept")
        performHC <- performHC && has_interceptX && (has_interceptX == has_interceptZ)
        indBetas <- if (performHC) {
            terms.labs_X <- attr(TermsX, "term.labels")
            terms.labs_Z <- attr(TermsZ, "term.labels")
            # check for time-varying covariates
            timeTerms <- grep(timeVar, colnames(X), fixed = TRUE)
            which_td <- unname(which(apply(X, 2, check_td, id = id)))
            all_TDterms <- unique(c(timeTerms, which_td))
            baseline <- seq_len(ncol(X))[-all_TDterms]
            #factorsX <- attr(TermsX, "factors")           
            c(list(baseline), lapply(colnames(Z)[-1L], find_positions, 
                                     nams2 = colnames(X)))
        }
    } else {
        performHC <- FALSE
        if (param == "shared-betasRE") {
            warning("\nit seems that the random effects design matrix is not a subset of the ",
                    "fixed effects design matrix. Argument 'param' is set to 'shared-RE'.")
            param <- "shared-RE"
        }
    }
    if (check_names && param == "shared-betasRE") {
        indBetasRE <- match(colnames(Z), colnames(X))
    }
    # transormation functions
    if (is.null(transFun)) {
        transFun.value <- function (x, data) x
        transFun.extra <- function (x, data) x
    } else {
        if (is.function(transFun)) {
            transFun.value <- transFun.extra <- transFun
        }
        if (is.list(transFun)) {
            transFun.value <- transFun$value
            transFun.extra <- transFun$extra
        }
        if (any(!names(formals(transFun.value)) %in% c("x", "data")))
            stop("\nincorrect specification of 'transFun' arguments.")
        if (any(!names(formals(transFun.extra)) %in% c("x", "data")))
            stop("\nincorrect specification of 'transFun' arguments.")
    }
    # put functions in a list
    hasScale <- inherits(try(densLong(y = y.long[1L], eta.y = lmeObject$fitted[1L], 
                                      log = FALSE, data = data), 
                             silent = TRUE), "try-error")
    Funs <- list(transFun.value = transFun.value, transFun.extra = transFun.extra,
                 densRE = densRE, densLong = densLong, hasScale = hasScale)
    # control values
    con <- list(adapt = FALSE, n.iter = 20000L, n.burnin = 3000L, n.thin = 10L, 
                n.adapt = 3000L, keepRE = TRUE, n.batch = 100L, priorVar = 100, 
                performHC = performHC, robust_baseHaz = FALSE, rescale_Bs.gammas = TRUE,
                knots = NULL, ObsTimes.knots = TRUE, 
                lng.in.kn = if (baseHaz == "P-splines") 15L else 5L, ordSpline = 4L, 
                seed = 1L, diff = 2L, 
                GQsurv = if (!estimateWeightFun) "GaussKronrod" else "GaussLegendre", 
                GQsurv.k = if (!estimateWeightFun) 15L else 17L,
                priorShapes = list(shape1 = dnorm, shape2 = dgamma, shape3 = dunif),
                verbose = TRUE, verbose2 = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (!any(namc == "n.thin"))
        con$n.thin <- max(1, floor(con$n.iter / 2000))
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    # construct desing matrices for longitudinal part for the hazard function
    data.id <- data[!duplicated(id), ]
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    if (param %in% c("td-value", "td-both")) {
        mfX.id <- model.frame(TermsX, data = data.id)
        mfZ.id <- model.frame(TermsZ, data = data.id)
        Xtime <- model.matrix(formYx, mfX.id)
        Ztime <- model.matrix(formYz, mfZ.id)
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame(terms(extraForm$fixed), data = data)
        TermsX.extra <- attr(mfX.extra, "terms")
        mfZ.extra <- model.frame(terms(extraForm$random), data = data)
        TermsZ.extra <- attr(mfZ.extra, "terms")
        mfX.extra.id <- model.frame(TermsX.extra, data = data.id)
        mfZ.extra.id <- model.frame(TermsZ.extra, data = data.id)
        Xtime.extra <- model.matrix(extraForm$fixed, mfX.extra.id)
        Ztime.extra <- model.matrix(extraForm$random, mfZ.extra.id)
        Xextra <- model.matrix(extraForm$fixed, mfX.extra)
        Zextra <- model.matrix(extraForm$random, mfZ.extra)
        long.extra <- as.vector(c(Xextra %*% fixef(lmeObject)[extraForm$indFixed]) +
            if (all(extraForm$indRandom > 0)) {
                rowSums(Zextra * b[id, extraForm$indRandom, drop = FALSE])
            } else {
                rep(0, nrow(Zextra))
            })
    }
    if (param == "td-value" || param == "shared-betasRE" || param == "shared-RE")
        long.extra <- NULL
    if (param == "td-extra")
        long <- NULL
    if (performHC) {
        data.idHC <- data.id
        data.idHC[[timeVar]] <- 1
        mfHC <- model.frame(TermsX, data = data.idHC)
        which.timeVar <- grep(timeVar, names(mfHC), fixed = TRUE)
        mfHC[which.timeVar] <- lapply(mfHC[which.timeVar], function (x) { x[] <- 1; x })
        XXtime <- model.matrix(formYx, mfHC)
    }
    # response vectors and design matrices
    y <- list(y = y.long, Time = Time, event = event, lag = lag, df.RE = df.RE, id = id,
              indBetas = if (check_names) indBetas, 
              indBetasRE = if (param == "shared-betasRE") indBetasRE, 
              LongFormat = LongFormat, offset = offset)
    if (typeSurvInf == "counting")
        y <- c(y, list(TimeL = TimeL, TimeR = TimeR, eventLong = eventLong, idT = idT,
                       typeSurvInf = typeSurvInf, anyLeftTrunc = anyLeftTrunc))
    x <- list(X = X, Z = Z, W = W)
    if (typeSurvInf == "counting") {
        wind <- tapply(idT, idT, function (x) rep(c(FALSE, TRUE), c(length(x) - 1, 1)))
        x$W <- x$W[unlist(wind, use.names = FALSE), , drop = FALSE]
    }
    x <- switch(param,
                "td-value" = c(x, list(Xtime = Xtime, Ztime = Ztime)),
                "td-extra" = c(x, list(Xtime.extra = Xtime.extra, Ztime.extra = Ztime.extra)),
                "td-both" = c(x, list(Xtime = Xtime, Ztime = Ztime,
                                      Xtime.extra = Xtime.extra, Ztime.extra = Ztime.extra)),
                "shared-betasRE" =, "shared-RE" =  x)
    if (performHC) {
        x <- c(x, list(XXtime = XXtime))
    }
    # construct desing matrices for longitudinal part for the survival function
    GQsurv <- if (con$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(con$GQsurv.k)
    wk <- GQsurv$wk
    sk <- GQsurv$sk
    K <- length(sk)
    P <- if (typeSurvInf == "counting" && anyLeftTrunc) (Time - TimeL) / 2 else Time / 2
    st <- if (typeSurvInf == "counting" && anyLeftTrunc) {
        outer(P, sk) + c(Time + TimeL) / 2
    } else {
        outer(P, sk + 1)
    }
    id.GK <- rep(seq_along(Time), each = K)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
    y <- c(y, list(id.GK = id.GK))
    x <- c(x, list(P = P, st = st, wk = wk))
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        Xs <- model.matrix(formYx, mfX)
        Zs <- model.matrix(formYz, mfZ)
        x <- c(x, list(Xs = Xs, Zs = Zs))
        if (estimateWeightFun) {
            P2 <- c(t(st)) / 2
            st2 <- outer(P2, sk + 1)
            id.GK2 <- rep(seq_len(nrow(data.id2)), each = K)
            data.id3 <- data.id2[id.GK2, ]
            data.id3[[timeVar]] <- pmax(c(t(st2)) - lag, 0)
            mfX <- model.frame(TermsX, data = data.id3)
            mfZ <- model.frame(TermsZ, data = data.id3)
            Xu <- model.matrix(formYx, mfX)
            Zu <- model.matrix(formYz, mfZ)
            x <- c(x, list(Xu = Xu, Zu = Zu, P2 = P2, st2 = st2))
            y <- c(y, list(id.GK2 = id.GK2))
            if (is.null(weightFun) || !is.function(weightFun)) {
                weightFun <- function (u, parms, t.max) {
                    num <- dnorm(x = u, mean = parms[1L], sd = parms[2L])
                    den <- pnorm(q = c(0, t.max), mean = parms[1L], sd = parms[2L])
                    num / (den[2L] - den[1L])
                }
            }
            Funs <- c(Funs, list(weightFun = weightFun))
        }
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame(TermsX.extra, data = data.id2)
        mfZ.extra <- model.frame(TermsZ.extra, data = data.id2)
        Xs.extra <- model.matrix(extraForm$fixed, mfX.extra)
        Zs.extra <- model.matrix(extraForm$random, mfZ.extra)
        x <- c(x, list(Xs.extra = Xs.extra, Zs.extra = Zs.extra))
    }
    # extra design matrices for the log approximated baseline hazard
    kn <- if (is.null(con$knots)) {
        if (baseHaz == "P-splines") {
            tt <- if (con$ObsTimes.knots) Time else Time[event == 1]
            pp <- quantile(tt, c(0.05, 0.95), names = FALSE)
            tail(head(seq(pp[1L], pp[2L], length.out = con$lng.in.kn), -1), -1)
        } else {
            pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
            pp <- tail(head(pp, -1), -1)
            tt <- if (con$ObsTimes.knots) Time else Time[event == 1]
            quantile(tt, pp, names = FALSE)            
        }
    } else {
        con$knots
    }
    kn <- kn[kn < max(Time)]
    rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
    con$knots <- rr
    W2 <- splineDesign(rr, Time, ord = con$ordSpline)
    if (any(colSums(W2) == 0))
        stop("\nsome of the knots of the B-splines basis are set outside the range",
             "\nof the observed event times.")
    W2s <- splineDesign(rr, c(t(st)), ord = con$ordSpline)
    if (typeSurvInf == "counting" && LongFormat) {
        TDind <- mapply(findInterval, x = split(st, row(st)), vec = split(TimeLl, idT), 
                        SIMPLIFY = FALSE)
        Ws <- do.call(rbind, mapply(function (x, i) x[i, , drop = FALSE], i = TDind,
                                    x = lapply(split(W, idT), matrix, ncol = ncol(W)), 
                                    SIMPLIFY = FALSE))
    } else {
        Ws <- NULL
    }
    x <- c(x, list(Ws = Ws, W2 = W2, W2s = W2s))
    # All data
    Data <- list(data = data, data.id = data.id, data.s = data.id2, 
                 data.u = if (estimateWeightFun) data.id3)
    # extract initial values
    D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
                 lmeObject$sigma^2)[[1]]
    sigma2 <- lmeObject$sigma^2
    initial.values <- list(betas = fixef(lmeObject), tau = 1/sigma2, b = b,
                           D = D, invD = invD)
    initSurv <- initSurvival(Time, event, id, W2, W2s, P, wk, id.GK, times, 
                             b, initial.values$betas, y$indBetas, 
                             x$W, baseHaz, con$diff, Data, param, 
                             if (param %in% c("td-value", "td-both")) long else NULL, 
                             long.extra, transFun.value, transFun.extra, 
                             vl = if (param %in% c("td-value", "td-both")) 
                                 transFun.value(c(Xtime %*% fixef(lmeObject)) + 
                                                    rowSums(Ztime * b), Data$data.id),
                             vls = if (param %in% c("td-value", "td-both")) 
                                 transFun.value(c(Xs %*% fixef(lmeObject)) + 
                                                    rowSums(Zs * b[id.GK, , drop = FALSE]), Data$data.s),
                             ex = if (param %in% c("td-extra", "td-both")) {
                                 iF <- extraForm$indFixed; iR <- extraForm$indRandom
                                 transFun.extra(c(Xtime.extra %*% fixef(lmeObject)[iF]) + 
                                                    rowSums(Ztime.extra * b[, iR, drop = FALSE]), Data$data.id)
                             },
                             exs = if (param %in% c("td-extra", "td-both")) {
                                 iF <- extraForm$indFixed; iR <- extraForm$indRandom
                                 transFun.extra(c(Xs.extra %*% fixef(lmeObject)[iF]) + 
                                                    rowSums(Zs.extra * b[id.GK, iR, drop = FALSE]), Data$data.s)
                             })
    initial.values$gammas <- initSurv$gammas
    initial.values$Bs.gammas <- initSurv$Bs.gammas
    initial.values$tauBs <- initSurv$tauBs
    initial.values$deltaBs <- 0.005
    initial.values$alphas <- initSurv$alphas
    initial.values$Dalphas <- initSurv$Dalphas
    if (estimateWeightFun) {
        w1 <- try(weightFun(0.1, c(0.1), 1), TRUE)
        w2 <- try(weightFun(0.1, c(0.1, 0.1), 1), TRUE)
        w3 <- try(weightFun(0.1, c(0.1, 0.1, 0.1), 1), TRUE)
        ind1 <- inherits(w1, "try-error") || is.na(w1)
        ind2 <- inherits(w2, "try-error") || is.na(w2)
        ind3 <- inherits(w3, "try-error") || is.na(w3)
        nshapes <- if (ind1 && ind2) {
            3
        } else if (ind2 || (!ind1 && w1 == w2)) {
            1
        } else if (ind3 || (!ind2 && w2 == w3)) {
            2
        } else 3
        initial.values$shapes <- rep(0.1, nshapes)
    }
    initial.values <- initial.values[!sapply(initial.values, is.null)]
    if (!is.null(init)) {
        lngths <- lapply(initial.values[(nam.init <- names(init))], length)
        if (!is.list(init) || !isTRUE(all.equal(lngths, lapply(init, length)))) {
            warning("'init' is not a list with elements numeric vectors of appropriate ",
                    "length; default starting values are used instead.\n")
        } else {
            initial.values[nam.init] <- init
            if (!is.matrix(initial.values$b))
                dim(initial.values$b) <- dim(b)
            if (!is.matrix(initial.values$D))
                dim(initial.values$D) <- dim(initial.values$invD) <- dim(D)            
        }
    }
    # default priors
    prs <- list(priorMean.betas = numeric(ncol(X)), 
                priorTau.betas = diag(1 / con$priorVar, ncol(X)),
                priorA.tau = (1/sigma2)^2 / 10, priorB.tau = (1/sigma2) / 10,
                priorR.invD = ncol(Z) * invD, priorK.invD = ncol(Z),
                priorMean.Bs.gammas = numeric(ncol(W2)), 
                priorTau.Bs.gammas = diag(10 / con$priorVar, ncol(W2)))
    if (baseHaz == "P-splines") {
        DD <- diag(ncol(W2))
        prs$priorTau.Bs.gammas <- crossprod(diff(DD, differences = con$diff)) + 1e-06 * DD
        #prs$priorA.tauBs <- 1
        #prs$priorB.tauBs <- 0.005
        prs$priorA.tauBs <- 1
        prs$priorB.tauBs <- if (con$robust_baseHaz) 1 else 0.005
        prs$priorA.deltaBs <- 1e-02
        prs$priorB.deltaBs <- 1e-02
    }
    if (!is.null(W)) {
        prs$priorMean.gammas <- numeric(ncol(W))
        prs$priorTau.gammas <- drop(diag(1 / con$priorVar, ncol(W)))
    }
    if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) {
        prs$priorMean.alphas <- numeric(length(initial.values$alphas))
        prs$priorTau.alphas <- drop(diag(10 / con$priorVar, length(initial.values$alphas)))
    }
    if (param %in% c("td-extra", "td-both")) {
        prs$priorMean.Dalphas <- numeric(length(initial.values$Dalphas))
        prs$priorTau.Dalphas <- drop(diag(10 / con$priorVar, length(initial.values$Dalphas)))
    }
    if (estimateWeightFun) {
        maxT <- max(Time) * 0.7
        prs$priorshape1 <- c(0, 5)
        prs$priorshape2 <- c(0.1, 0.1)
        prs$priorshape3 <- c(-maxT, maxT)
    }
    if (!is.null(priors)) {
        lngths <- lapply(prs[(nam.prs <- names(priors))], length)
        if (!is.list(priors) || !isTRUE(all.equal(lngths, lapply(priors, length)))) {
            warning("'priors' is not a list with elements numeric vectors of appropriate ",
                    "length; default priors are used instead.\n")
        } else {
            prs[nam.prs] <- priors
        }
    }
    # covariance matrices of parameters from separate fits
    Covs <- list(betas = vcov(lmeObject), b = Cov.postRE, 
                 gammas = if (!is.null(W)) initSurv$cov.gammas, Bs.gammas = initSurv$cov.Bs.gammas,
                 alphas = initSurv$cov.alphas, Dalphas = initSurv$cov.Dalphas)
    # delete unused objects to free memory
    namsDelObjs <- c(names(x), "y.long", "Time", "event", "data.id2", "b", "D", "invD", "Cov.postRE",
                     "id", "id.GK", "kn", "i", "K", "Z.i")
    if (param %in% c("td-value", "td-both")) {
        namsDelObjs <- c(namsDelObjs, "mfX", "mfZ", "mfX.id", "mfZ.id")
    }
    if (param %in% c("td-extra", "td-both")) {
        namsDelObjs <- c(namsDelObjs, "mfX.extra", "mfZ.extra", "mfX.extra.id", "mfZ.extra.id",
                         "Xextra", "Zextra")
    }    
    rm(list = namsDelObjs); gc()
    # joint model fit
    model <- MCMCfit(y, x, param, extraForm, baseHaz, estimateWeightFun, initial.values, prs, 
                     scales, Funs, Covs, Data, con, df.RE)
    names.betas <- names(fixef(lmeObject))
    names.D <- colnames(getVarCov(lmeObject))
    names.gammas <- names(coef(survObject))
    names.Bs.gammas <- paste0("Bs.gammas", seq_along(model$postMeans$Bs.gammas))
    names.alphas <- if (param %in% c("td-value", "td-both")) {
        if ((na <- length(model$postMeans$alphas)) == 1) "Assoct" else {
            nm <- colnames(transFun.value(1, data.id))
            if (all(nm[-1] == "")) paste0("Assoct", seq_len(na)) else c("Assoct", paste0("Assoct:", nm[-1]))
        }
    } else if (param %in% c("shared-betasRE", "shared-RE")) {
        paste0("Assoct:", names.D)
    }
    names.Dalphas <- if (param %in% c("td-extra", "td-both")) {
        if ((nda <- length(model$postMeans$Dalphas)) == 1) "AssoctE" else {
            nm <- colnames(transFun.extra(1, data.id))
            if (all(nm[-1] == "")) paste0("AssoctE", seq_len(nda)) else c("AssoctE", paste0("AssoctE:", nm[-1]))
        }
    }
    model <- fixNames(model, names.betas, names.D, names.gammas, names.Bs.gammas, 
                      names.alphas, names.Dalphas, 
                      names.shapes = if (estimateWeightFun) 
                          paste0("shape", seq_along(model$postMeans$shapes)),
                      names.id = row.names(ranef(lmeObject)))
    out <- c(model, list(x = x, y = y, Data = Data, Funs = Funs,
                         Terms = list(termsYx = TermsX, termsYz = TermsZ, 
                                      termsT = survObject$terms,
                                      termsYx.extra = if (param %in% c("td-extra", "td-both")) TermsX.extra, 
                                      termsYz.extra = if (param %in% c("td-extra", "td-both")) TermsZ.extra),
                         Forms = list(formYx = formYx, formYz = formYz, 
                                      formT = formT, extraForm = extraForm),
                         timeVar = timeVar, control = con, densLongCheck = densLongCheck,
                         param = param, priors = prs, baseHaz = baseHaz, df.RE = df.RE, 
                         estimateWeightFun = estimateWeightFun,
                         call = cl))
    class(out) <- "JMbayes"
    out
}
