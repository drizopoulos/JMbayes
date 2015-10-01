predict.JMbayes <- function (object, newdata, type = c("Marginal", "Subject"),
    interval = c("none", "confidence", "prediction"), level = 0.95, idVar = "id", 
    FtTimes = NULL, last.time = NULL, LeftTrunc_var = NULL, M = 300, returnData = FALSE, 
    scale = 1.6, weight = rep(1, nrow(newdata)), invlink = NULL, seed = 1, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    type <- match.arg(type)
    interval <- match.arg(interval)
    if (type == "Marginal") {
        TermsX <- delete.response(object$Terms$termsYx)
        mf <- model.frame(TermsX, data = newdata)
        form <- reformulate(attr(TermsX, "term.labels"))
        X <- model.matrix(form, data = mf)
        out <- c(X %*% object$postMeans$betas)
        if (!is.null(invlink))
            out <- invlink(out)
        names(out) <- row.names(newdata)
        if (interval == "prediction") {
            warning("\nfor type = 'Marginal' only confidence intervals are calculated.")
            interval <- "confidence"
        }
        if (interval == "confidence") {
            betas <- object$mcmc$betas
            preds <- X %*% t(betas)
            se.fit <- apply(preds, 1, sd)
            alpha <- 1 - level
            low <- apply(preds, 1, quantile, probs = alpha/2)
            up <- apply(preds, 1, quantile, probs = 1 - alpha/2)
            names(se.fit) <- names(low) <- names(up) <- row.names(newdata)
            out <- list(pred = out, se.fit = se.fit, low = low, upp = up)
            if (!is.null(invlink)) {
                out$low <- invlink(out$low)
                out$upp <- invlink(out$upp)
            }
        }
        if (returnData) {
            out <- if (is.list(out)) 
                cbind(newdata, do.call(cbind, out))
            else
                cbind(newdata, pred = out)
        }
    } else {
        if (interval == "prediction" && !is.null(invlink)) {
            interval <- "confidence"
            warning("\nprediction itervals are currently not calculated when 'invlink' is specified.\n")
        }
        if (is.null(newdata[[idVar]]))
            stop("'idVar' not in 'newdata.\n'")
        timeVar <- object$timeVar
        df.RE <- object$y$df.RE
        param <- object$param
        densLong <- object$Funs$densLong
        hasScale <- object$Funs$hasScale
        anyLeftTrunc <- object$y$anyLeftTrunc
        densRE <- object$Funs$densRE
        transFun.value <- object$Funs$transFun.value
        transFun.extra <- object$Funs$transFun.extra
        extraForm <- object$Forms$extraForm
        indFixed <- extraForm$indFixed
        indRandom <- extraForm$indRandom
        indBetas <- object$y$indBetas
        TermsX <- object$Terms$termsYx
        TermsZ <- object$Terms$termsYz
        TermsX.extra <- object$Terms$termsYx.extra
        TermsZ.extra <- object$Terms$termsYz.extra
        mfX <- model.frame(TermsX, data = newdata)
        mfZ <- model.frame(TermsZ, data = newdata)
        formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
        formYz <- object$Forms$formYz
        estimateWeightFun <- object$estimateWeightFun
        weightFun <- object$Funs$weightFun
        max.time <- max(object$y$Time)
        na.ind <- as.vector(attr(mfX, "na.action"))
        na.ind <- if (is.null(na.ind)) {
            rep(TRUE, nrow(newdata))
        } else {
            !seq_len(nrow(newdata)) %in% na.ind
        }
        id <- as.numeric(unclass(newdata[[idVar]]))
        id <- id. <- match(id, unique(id))
        id <- id[na.ind]
        y <- model.response(mfX)
        X <- model.matrix(formYx, mfX)
        Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
        TermsT <- object$Terms$termsT
        data.id <- newdata[tapply(row.names(newdata), id, tail, n = 1L), ]
        data.s <- data.id[rep(1:nrow(data.id), each = 15L), ]
        idT <- data.id[[idVar]]
        idT <- match(idT, unique(idT))
        ids <- data.s[[idVar]]
        ids <- match(ids, unique(ids))
        mfT <- model.frame(delete.response(TermsT), data = data.id)        
        formT <- if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
            tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
            reformulate(attr(tt, "term.labels"))
        } else {
            tt <- attr(delete.response(TermsT), "term.labels")
            if (length(tt)) reformulate(tt) else reformulate("1")
        }
        W <- model.matrix(formT, mfT)[, -1, drop = FALSE]
        obs.times <- split(newdata[[timeVar]], id.)
        last.time <- if (is.null(last.time)) {
            tapply(newdata[[timeVar]], id., tail, n = 1L)
        } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
            last.time
        } else {
            stop("\nnot appropriate value for 'last.time' argument.")
        }
        times <- object$Data$data[[timeVar]]
        times.to.pred <- if (is.null(FtTimes)) {
            lapply(last.time, 
                function (t) seq(t, max(times) + 0.1 * mad(times), length = 25L))
        } else {
            if (!is.list(FtTimes) || length(FtTimes) != length(last.time))
                rep(list(FtTimes), length(last.time))
            else
                FtTimes
        }
        TimeL <- if (!is.null(anyLeftTrunc) && anyLeftTrunc) {
            if (is.null(LeftTrunc_var) || is.null(newdata[[LeftTrunc_var]])) {
                warning("The original joint model was fitted in a data set with left-",
                        "truncation and\nargument 'LeftTrunc_var' of predict.JMbayes() has not ", 
                        "been specified.\n")
            }
            TimeL <- newdata[[LeftTrunc_var]]
            tapply(TimeL, id, head, n = 1)
        }
        n <- length(object$y$Time)
        n.tp <- length(last.time)
        ncx <- ncol(X)
        ncz <- ncol(Z)
        ncww <- ncol(W)
        if (ncww == 0)
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
        list.thetas <- list(betas = betas, sigma = sigma, gammas = gammas, 
                            alphas = alphas, Dalphas = Dalphas, shapes = shapes, 
                            Bs.gammas = Bs.gammas, D = D)
        list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
        thetas <- unlist(as.relistable(list.thetas))
        Var.thetas <- vcov(object)
        environment(log.posterior.b) <- environment(ModelMats) <- environment()
        # construct model matrices to calculate the survival functions
        obs.times.surv <- split(data.id[[timeVar]], idT)
        survMats.last <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            survMats.last[[i]] <- ModelMats(last.time[i], ii = i)
        }
        data.id2 <- newdata[!duplicated(id), ]
        data.id2 <- data.id2[rep(1:nrow(data.id2), 
            sapply(times.to.pred, length)), ]
        data.id2[[timeVar]] <- unlist(times.to.pred)
        mfXpred <- model.frame(TermsX, data = data.id2)
        mfZpred <- model.frame(TermsZ, data = data.id2)
        Xpred <- model.matrix(formYx, mfXpred)
        Zpred <- model.matrix(formYz, mfZpred)
        id2 <- as.numeric(unclass(data.id2[[idVar]]))
        id2 <- match(id2, unique(id2))
        # calculate the Empirical Bayes estimates and their (scaled) variance
        modes.b <- matrix(0, n.tp, ncz)
        invVars.b <- Vars.b <- vector("list", n.tp)
        set.seed(seed)
        for (i in seq_len(n.tp)) {
            betas.new <- betas
            sigma.new <- sigma
            D.new <- D
            gammas.new <- gammas
            alphas.new <- alphas
            Dalphas.new <- Dalphas
            shapes.new <- shapes
            Bs.gammas.new <- Bs.gammas
            ff <- function (b, y, tt, mm, i) 
                -log.posterior.b(b, y, Mats = tt, ii = i)
            opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, i = i, 
                method = "BFGS", hessian = TRUE), TRUE)
            if (inherits(opt, "try-error")) {
                gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, i = i)
                opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last,
                    i = i, method = "BFGS", hessian = TRUE)
            }
            modes.b[i, ] <- opt$par
            invVars.b[[i]] <- opt$hessian/scale
            Vars.b[[i]] <- scale * solve(opt$hessian)        
        }
        res <- vector("list", M)
        success.rate <- matrix(FALSE, M, n.tp)
        b.old <- b.new <- modes.b
        if (n.tp == 1)
            dim(b.old) <- dim(b.new) <- c(1L, ncz)    
        mcmc <- object$mcmc
        mcmc <- mcmc[names(mcmc) != "b"]
        if (M > nrow(mcmc$betas)) {
            warning("'M' cannot be set greater than ", nrow(mcmc$betas))
            M <- nrow(mcmc$betas)
            out <- vector("list", M)
            success.rate <- matrix(FALSE, M, n.tp)
        }
        samples <- sample(nrow(mcmc$betas), M)
        mcmc[] <- lapply(mcmc, function (x) x[samples, , drop = FALSE])
        proposed.b <- mapply(rmvt, mu = split(modes.b, row(modes.b)), Sigma = Vars.b, 
                             MoreArgs = list(n = M, df = 4), SIMPLIFY = FALSE)
        dmvt.proposed <- mapply(dmvt, x = proposed.b, mu = split(modes.b, row(modes.b)),
                                Sigma = Vars.b, MoreArgs = list(df = 4, log = TRUE), SIMPLIFY = FALSE)
        for (m in 1:M) {
            # Step 1: simulate new parameter values
            betas.new <- mcmc$betas[m, ]
            if (hasScale)
                sigma.new <- mcmc$sigma[m, ]
            if (!is.null(W))
                gammas.new <- mcmc$gammas[m, ]
            if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
                alphas.new <- mcmc$alphas[m, ]
            if (param %in% c("td-extra", "td-both"))
                Dalphas.new <- mcmc$Dalphas[m, ]
            if (estimateWeightFun)
                shapes.new <- mcmc$shapes[m, ]
            D.new <- mcmc$D[m, ]; dim(D.new) <- dim(D)
            Bs.gammas.new <- mcmc$Bs.gammas[m, ]
            y.new <- vector("list", n.tp)
            for (i in seq_len(n.tp)) {
                # Step 2: simulate new random effects values
                p.b <- proposed.b[[i]][m, ]
                dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], invSigma = invVars.b[[i]], df = 4, log = TRUE)
                dmvt.prop <- dmvt.proposed[[i]][m]
                a <- min(exp(log.posterior.b(p.b, y, survMats.last, ii = i) + dmvt.old - 
                        log.posterior.b(b.old[i, ], y, survMats.last, ii = i) - dmvt.prop), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (!is.na(ind) && ind)
                    b.new[i, ] <- p.b
                # Step 3: compute future Ys
                Xpred.i <- Xpred[id2 == i, , drop = FALSE]
                Zpred.i <- Zpred[id2 == i, , drop = FALSE]                
                mu.i <- as.vector(c(Xpred.i %*% betas.new) + 
                    rowSums(Zpred.i * rep(b.new[i, ], each = nrow(Zpred.i))))
                if (!is.null(invlink))
                    mu.i <- invlink(mu.i)
                y.new[[i]] <- if (interval == "confidence") weight[i] * mu.i else 
                    if (interval == "prediction") {
                            weight[i] * rnorm(length(mu.i), mu.i, sigma.new)
                    }
            }
            b.old <- b.new
            res[[m]] <- y.new
        }
        oo <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            oo[[i]] <- do.call(rbind, sapply(res, "[", i))
        }
        out <- as.vector(c(Xpred %*% betas) + rowSums(Zpred * modes.b[id2, , drop = FALSE]))
        if (!is.null(invlink))
            out <- invlink(out)
        out <- weight[i] * out
        if (interval %in% c("confidence", "prediction")) {
            alpha <- 1 - level
            se.fit <- lapply(oo, function (m) apply(as.matrix(m), 2, sd))
            f1 <- function (mat) apply(mat, 2, quantile, probs = alpha/2)
            f2 <- function (mat) apply(mat, 2, quantile, probs = 1 - alpha/2)
            low <- lapply(oo, f1) 
            up <- lapply(oo, f2)
            out <- list(pred = out, se.fit = unlist(se.fit), 
                low = unlist(low), upp = unlist(up), all.vals = oo)
            if (!is.null(invlink)) {
                out$low <- invlink(out$low)
                out$upp <- invlink(out$upp)
            }
        }
        if (returnData) {
            newdata$pred <- c(X %*% betas) + rowSums(Z * modes.b[id, ])
            out <- if (is.list(out)) {
                newdata$upp <- newdata$low <- newdata$se.fit <- NA
                rbind(newdata, cbind(data.id2, do.call(cbind, out[-5])))
            } else {
                rbind(newdata, cbind(data.id2, pred = out))
            }
        } else
            attr(out, "time.to.pred") <- times.to.pred
    }
    rm(list = ".Random.seed", envir = globalenv())
    class(out) <- c(class(out), "predict.JMbayes")
    out
}
