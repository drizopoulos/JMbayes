mvglmer <- function (formulas, data, families, engine = c("JAGS", "STAN"), 
                     overdispersion = FALSE, priors = NULL, init = NULL, 
                     control = NULL, ...) {
    cl <- match.call()
    engine <- match.arg(engine)
    if (!is.list(families))
        stop("'families' must be a list of family objects.")
    # depending on the input of the user, set families to the corresponding
    # base R functions
    families[] <- lapply(families, function (family) {
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family))
            family <- family()
        family
    })
    # using the formulas and the data extract and constuct the objects required to
    # pass to JAGS
    components <- lapply(unname(formulas), extractFrames, data = data)
    # perform some checks
    if (!all((ns <- sapply(components, `[[`, 'n')) == components[[1L]][['n']])) {
        stop("it seems that the number of subjects differ between the responses. More ",
             "specifically, according to the data the number of subjects per response is ",
             paste(paste(sapply(components, `[[`, 'respVar'), ns, sep = " = "),
                   collapse = ", "), ".\n")
    }
    components <- unlist(components, recursive = FALSE)
    n_outcomes <- length(formulas)
    names(components) <- paste0(names(components),
                                rep(seq_len(n_outcomes),
                                    each = length(components) / n_outcomes))
    colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
    colmns_nHC <- components[grep("colmns_nHC", names(components), fixed = TRUE)]
    seq_outcomes <- seq_len(n_outcomes)
    nams_vars <- c("N", "id", "Z", "Xhc", "ncx", "y")
    vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
    if (any(ind_td <- sapply(colmns_nHC, length))) {
        vars <- c(vars, paste0("X", which(ind_td > 0)))
    }
    Data <- c(list(n = components$n1), components[vars])
    Data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
    RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                      sq = components[grep("ncz", names(components), fixed = TRUE)],
                      incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                      SIMPLIFY = FALSE)
    names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
    Data <- c(Data, RE_inds, unlist(colmns_HC, recursive = FALSE), colmns_nHC)
    # control
    con <- list(n.processors = parallel::detectCores() - 1, n.chains = 2,
                working.directory = getwd(), clear.model = TRUE,
                seed = 1L, optimize_only = FALSE, verbose = FALSE)
    if (engine == "JAGS") {
        con$n.iter <- 28000L
        con$n.burnin <- 3000L
        con$n.thin <- 50L
        con$n.adapt <- 3000L 
    } else {
        con$n.iter <- 1000
        con$n.warmup <- floor(con$n.iter / 2)
        con$n.thin <- 1
        con$adapt_delta <- 0.8
    }
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (!any(namc == "n.thin")) {
        con$n.thin <- if (engine == "JAGS") {
            max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
        } else {
            max(1, floor((con$n.iter - con$n.warmup) * con$n.chains / 1000))
        }
    }
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ######################################################################################
    # Priors
    if (engine == "JAGS") {
        prs <- list(priorR_D = diag(rep(as.numeric(NA), Data$n_RE), Data$n_RE),
                    priorK_D = Data$n_RE + 1, A_RD = 0.5, B_RD = 0.01,
                    tau_half_cauchy = 0.1)
        pr_taus_betas <- rep(list(0.01), n_outcomes)
        names(pr_taus_betas) <- paste0("tau_betas", seq_len(n_outcomes))
        prs <- c(prs, pr_taus_betas)
        if (any(sapply(families, function (x) x$family == "gaussian")) || overdispersion) {
            prs$A_tau <- 0.01
            prs$B_tau <- 0.01
        }
    } else {
        prs <- list(scale_sigmas = 5, scale_diag_D = 3, lkj_shape = 2,
                    priorK_D = Data$n_RE + 1)
        pr_scale_betas <- rep(list(10), n_outcomes)
        names(pr_scale_betas) <- paste0("scale_betas", seq_len(n_outcomes))
        prs <- c(prs, pr_scale_betas)
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
    Data <- c(Data, prs)
    ######################################################################################
    # write model
    model_name <- paste0("mv_glmer", sample(1e06, 1), if (engine == "JAGS") ".txt" else ".stan")
    if (engine == "JAGS") {
        cat(build_model(families, seq_along(families), colmns_HC, colmns_nHC, overdispersion,
                        Data$n_RE), file = file.path(con$working.directory, model_name))
    } else {
        cat(data_part(families, lapply(colmns_HC, length), lapply(colmns_nHC, length), 
                      Data$n_RE, colmns_HC, colmns_nHC),
            parameters(families, Data$n_RE),
            transformed_parameters(families, colmns_HC, colmns_nHC, RE_inds),
            model(families, Data$n_RE), generated_quantities(Data$n_RE),
            file = file.path(con$working.directory, model_name))
    }
    #closeAllConnections()
    # parameters to save
    params <- paste0('betas', seq_len(n_outcomes))
    if (any(ind_gs <- sapply(families, function (x) x$family == "gaussian"))) {
        params <- c(params, paste0("sigma", which(ind_gs)))
    }
    params <- if (engine == "JAGS") c(params, "inv_D", "b") else c(params, "D", "b")
    inits <- function () {
        ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
                        rnorm, sd = 0.1)
        names(ints) <- paste0('betas', seq_len(n_outcomes))
        ints$u <- drop(matrix(rnorm(Data$n * Data$n_RE), Data$n,
                              Data$n_RE))
        ints$inv_D <- ints$D <- if (Data$n_RE > 1) diag(Data$n_RE) else 1
        if (any(ind_gs)) {
            nms <- which(ind_gs)
            taus <- rep(list(1), length(nms))
            names(taus) <- paste0("tau", nms)
            ints <- c(ints, taus)
        }
        ints
    }
    fit <- if (engine == "JAGS") {
        jagsUI::jags(data = Data, inits = inits, parameters.to.save = params,
                     model.file = file.path(con$working.directory, model_name),
                     parallel = con$n.processors > 1, n.chains = con$n.chains,
                     n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                     n.thin = con$n.thin, seed = con$seed, verbose = con$verbose)
    } else {
        options(mc.cores = con$n.chains)
        if (con$optimize_only) {
            model <- rstan::stan_model(file = file.path(con$working.directory, model_name))
            if (is.null(init)) {
                init <- "random"
            } else {
                if (!is.list(init)) {
                    stop("'init' should be a list with appropriate names; run once ",
                         "mvglmer() and see the components of object$par.")
                }
            }
            out <- rstan::optimizing(model, data = Data, hessian = TRUE, as_vector = FALSE,
                                     seed = con$n.iter, init = init)
            if (con$clear.model) {
                file.remove(file.path(con$working.directory, model_name))
            }
            return(out)
        }
        out <- rstan::stan(file = file.path(con$working.directory, model_name), data = Data, 
                    pars = params, iter = con$n.iter, chains = con$n.chains, 
                    thin = con$n.thin, seed = con$seed, 
                    control = list('adapt_delta' = con$adapt_delta))
        sims.list <- lapply(lapply(params, extract, object = out, permuted = FALSE), bind_chains)
        sims.list[] <- lapply(sims.list, function (x) 
            if (length(dim(x)) == 1) as.matrix(x) else x)
        names(sims.list) <- params
        sims.list[['D']] <- fix_D(sims.list[['D']])
        sims.list[['b']] <- fix_b(sims.list[['b']], Data$n_RE)
        splts <- rep(seq_along(sims.list), sapply(sims.list, function (x) prod(dim(x)[-1])))
        rhats <- head(rstan::summary(out)$summary[, "Rhat"], -1)
        Rhat <- split(rhats, splts)
        names(Rhat) <- params
        list(sims.list = sims.list, Rhat = Rhat,
             mcmc.info = list(n.chains = con$n.chains, n.thin = con$n.thin,
                              n.warmup = con$n.warmup,
                              n.samples = nrow(sims.list[["betas1"]]),
                              elapsed.mins = sum(get_elapsed_time(out)) / 60), 
             DIC = NULL, pD = NULL)
    }
    if (con$clear.model) {
        file.remove(file.path(con$working.directory, model_name))
    }
    out <- list(mcmc = fit$sims.list, components = components, data = data,
                families = families, control = con, mcmc.info = fit$mcmc.info,
                DIC = fit$DIC, pD = fit$pD, Rhat = fit$Rhat,
                priors = prs, engine = engine)
    if (Data$n_RE == 1) {
        if (engine == "JAGS") 
            out$mcmc$inv_D <- array(out$mcmc$inv_D, c(length(out$mcmc$inv_D), 1, 1))
        else 
            out$mcmc$D <- array(out$mcmc$D, c(length(out$mcmc$D), 1, 1))
        out$mcmc$b <- array(out$mcmc$b, c(nrow(out$mcmc$b), ncol(out$mcmc$b), 1))
    }
    if (engine == "JAGS") {
        out$mcmc$D <- out$mcmc$inv_D
        for (i in seq_len(nrow(out$mcmc$betas1))) {
            out$mcmc$D[i, , ] <- solve(out$mcmc$D[i, , ])
        }
    } else {
        out$mcmc$inv_D <- out$mcmc$D
        for (i in seq_len(nrow(out$mcmc$betas1))) {
            out$mcmc$inv_D[i, , ] <- solve(out$mcmc$D[i, , ])
        }
    }
    # fix names
    #pat <- paste0("^X", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
    Xnams <- lapply(components[grep("^X[0-9]", names(components))], colnames)
    for (i in seq_along(Xnams)) {
        colnames(out$mcmc[[paste0("betas", i)]]) <- Xnams[[i]]
    }
    #pat <- paste0("^Z", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
    Znams <- lapply(components[grep("^Z[0-9]", names(components))], colnames)
    Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
    dimnames(out$mcmc$D) <- dimnames(out$mcmc$inv_D) <- list(NULL, Znams, Znams)
    dimnames(out$mcmc$b) <- list(NULL, NULL, Znams)
    # calculate statistics
    summary_fun <- function (FUN, ...) {
        out <- lapply(out$mcmc, function (x) {
            if (!is.null(dim(x)) && length(dim(x)) > 1) {
                d <- if (is.matrix(x)) 2L else c(2L, 3L)
                apply(x, d, FUN, ...)
            } else {
                FUN(x, ...)
            }
        })
        out[!sapply(out, is.null)]
    }
    out$postMeans <- summary_fun(mean, na.rm = TRUE)
    out$postModes <- summary_fun(modes)
    out$EffectiveSize <- summary_fun(effectiveSize)
    out$StErr <- summary_fun(stdErr)
    out$StDev <- summary_fun(sd, na.rm = TRUE)
    out$CIs <- summary_fun(quantile, probs = c(0.025, 0.975))
    out$Pvalues <- summary_fun(computeP)
    out$call <- cl
    class(out) <- "mvglmer"
    out
}

summary.mvglmer <- function (object, ...) {
    families <- object$families
    n_outcomes <- length(families)
    components <- object$components
    extract_components <- function (nam) {
        components[grep(nam, names(components), fixed = TRUE)]
    }
    respVars <- unlist(extract_components("respVar"), use.names = FALSE)
    descrpt <- data.frame(" " = unlist(extract_components("N"), use.names = FALSE),
                          row.names = respVars, check.rows = FALSE, check.names = FALSE)
    out <- list(n = components$n1, descrpt = descrpt, D = object$postMeans$D,
                families = families, respVars = respVars,
                control = object$control, mcmc.info = object$mcmc.info,
                DIC = object$DIC, pD = object$pD, call = object$call, engine = object$engine)
    for (i in seq_len(n_outcomes)) {
        out[[paste0("Outcome", i)]] <- data.frame("PostMean" = object$postMeans[[paste0("betas", i)]],
                                                  "StDev" = object$StDev[[paste0("betas", i)]],
                                                  "StErr"= object$StErr[[paste0("betas", i)]],
                                                  "2.5%" = object$CIs[[paste0("betas", i)]][1, ],
                                                  "97.5%" = object$CIs[[paste0("betas", i)]][2, ],
                                                  "P" = object$Pvalues[[paste0("betas", i)]],
                                                  "Rhat" = object$Rhat[[paste0("betas", i)]],
                                                  row.names = names(object$postMeans[[paste0("betas", i)]]),
                                                  check.names = FALSE)
        if (families[[i]][["family"]] == "gaussian") {
            D <- data.frame("PostMean" = object$postMeans[[paste0("sigma", i)]],
                            "StDev" = object$StDev[[paste0("sigma", i)]],
                            "StErr"= object$StErr[[paste0("sigma", i)]],
                            "2.5%" = object$CIs[[paste0("sigma", i)]][1],
                            "97.5%" = object$CIs[[paste0("sigma", i)]][2],
                            "P" = object$Pvalues[[paste0("sigma", i)]],
                            "Rhat" = object$Rhat[[paste0("sigma", i)]],
                            row.names = "sigma", check.names = FALSE)
            out[[paste0("Outcome", i)]] <- rbind(out[[paste0("Outcome", i)]], D)
        }
    }
    class(out) <- "summary.mvglmer"
    out
}

print.summary.mvglmer <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("Data Descriptives:")
    cat("\nNumber of Groups:", x$n)
    cat("\nNumber of Observations:\n")
    print(x$descrpt)
    cat("\n")
    if (!is.null(x$DIC)){
        model.sum <- data.frame(DIC = x$DIC, pD = x$pD, row.names = "")
        print(model.sum)
    }
    cat("\nRandom-effects covariance matrix:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- as.data.frame(round(rbind(sds), digits))
            names(dat) <- "StdDev"
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- rbind(mat)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- c(dimnames(D)[[1]])
        }
    } else {
        dat <- data.frame("StdDev" = c(sds, x$sigma),
                          row.names = if (!is.null(x$sigma)) c(rownames(D), "Residual") else rownames(D),
                          check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    n_outcomes <- length(x$families)
    for (i in seq_len(n_outcomes)) {
        cat("\nOutcome:", x$respVars[i],"\n")
        print(round(x[[paste0("Outcome", i)]], digits))
    }
    cat("\nMCMC summary:\n")
    tt <- x$mcmc.info$elapsed.mins
    cat("engine:", x$engine, 
        "\niterations:", x$control$n.iter, 
        if (x$engine == "JAGS") 
            paste("\nadapt:", x$control$n.adapt,
                  "\nburn-in:", x$control$n.burnin)
        else 
            paste("\nwarmup:", x$control$n.warmup), 
        "\nthinning:", x$control$n.thin,
        "\ntime:", if (tt > 60) round(tt/60, 1) else round(tt, 1),
        if (tt > 60) "hours" else "min")
    cat("\n\n")
    invisible(x)
}

plot.mvglmer <- function (x, which = c("trace", "autocorr", "density"),
                            param = c("betas", "sigma", "D"),
                            ask = TRUE, ...) {
    if (!inherits(x, "mvglmer"))
        stop("Use only with 'mvglmer' objects.\n")
    which <- match.arg(which)
    if (which %in% c("trace", "density", "autocorr")) {
        param <- match.arg(param, several.ok = TRUE)
        if (any(param == "D")) {
            keepD <- lower.tri(x$postMeans$D, TRUE)
            x$mcmc$D <- t(apply(x$mcmc$D, 1, c))[, c(keepD)]
            dnams <- which(keepD, arr.ind = TRUE)
            colnames(x$mcmc$D) <- paste0("D[", dnams[, 1], ", ", dnams[, 2], "]")
        }
        if (any(param == "tauBs")) {
            colnames(x$mcmc$tauBs) <- "tauBs"
        }
        which_parms <- unlist(sapply(param,
                                     function (pat) grep(paste0("^", pat), names(x$mcmc))),
                              use.names = FALSE)
        pp <- do.call(cbind, x$mcmc[which_parms])
        nams <- colnames(pp)
        op <- if (ask) par(mfrow = c(2, 2), ask = ask) else par(mfrow = c(4, 2))
        if (which == "trace") {
            for (i in 1:ncol(pp))
                plot(pp[, i], type = "l", xlab = "iterations", ylab = nams[i])
        } else if (which == "density") {
            for (i in 1:ncol(pp)) {
                bw <- bw.SJ(pp[, i]) * 1.5
                plot(density(pp[, i], bw = bw), xlab = nams[i],
                     main = paste("Density of", nams[i]))
            }
        } else {
            for (i in 1:ncol(pp))
                acf(pp[, i], ylab = nams[i], main = paste("Series", nams[i]))
        }
        par(op)
    }
    invisible()
}

fixef.mvglmer <- function (object, ...) {
    if (!inherits(object, "mvglmer"))
        stop("Use only with 'mvglmer' objects.\n")
    comps <- object$components
    nams_outcomes <- unlist(comps[grep("respVar", names(comps), fixed = TRUE)], 
                            use.names = FALSE)
    pMeans <- object$postMeans
    betas <- pMeans[grep("betas", names(pMeans), fixed = TRUE)]
    names(betas) <- nams_outcomes
    betas
}

bind_chains <- function (ar) {
    d <- dim(ar)
    e <- seq_len(d[2L]) * d[1L]
    s <- c(1, head(e, -1) + 1)
    ind <- mapply(seq, from = s, to = e, SIMPLIFY = FALSE)
    m <- array(0.0, c(d[1L] * d[2L], d[3L]))
    for (i in seq_len(d[2L])) {
        m[ind[[i]], ] <- ar[, i, ]
    }
    colnames(m) <- dimnames(ar)[[3]]
    m
}

fix_D <- function (D) {
    d <- dim(D)
    k <- round(sqrt(d[2L]))
    m <- array(0.0, c(d[1L], k, k))
    for (i in seq_len(d[1L]))
        m[i, , ] <- matrix(D[i, ], k, k)
    m
}

fix_b <- function (b, n_RE) {
    d <- dim(b)
    n <- round(d[2L] / n_RE)
    m <- array(0.0, c(d[1L], n, n_RE))
    for (i in seq_len(d[1L]))
        m[i, , ] <- matrix(b[i, ], n, n_RE)
    m
}
