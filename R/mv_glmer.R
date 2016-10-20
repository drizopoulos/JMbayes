mvglmer <- function (formulas, data, families, overdispersion = FALSE,
                      priors = NULL, control = NULL, ...) {
    cl <- match.call()
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
    jags_vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
    if (any(ind_td <- sapply(colmns_nHC, length))) {
        jags_vars <- c(jags_vars, paste0("X", which(ind_td > 0)))
    }
    JAGS_data <- c(list(n = components$n1), components[jags_vars])
    JAGS_data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
    RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                      sq = components[grep("ncz", names(components), fixed = TRUE)],
                      incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                      SIMPLIFY = FALSE)
    names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
    JAGS_data <- c(JAGS_data, RE_inds)
    # control
    con <- list(n.processors = detectCores() - 1, n.chains = 2,
                n.iter = 28000L, n.burnin = 3000L, n.thin = 50L,
                n.adapt = 3000L, working.directory = getwd(), clear.model = TRUE,
                seed = 1L, verbose = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (!any(namc == "n.thin"))
        con$n.thin <- max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ######################################################################################
    # Priors
    prs <- list(priorR.D = diag(rep(NA, JAGS_data$n_RE), JAGS_data$n_RE),
                priorK.D = JAGS_data$n_RE + 1, A_R.D = 0.5, B_R.D = 0.01,
                tau_half_cauchy = 0.1)
    pr_taus_betas <- rep(list(0.001), n_outcomes)
    names(pr_taus_betas) <- paste0("tau_betas", seq_len(n_outcomes))
    prs <- c(prs, pr_taus_betas)
    if (any(sapply(families, function (x) x$family == "gaussian")) || overdispersion) {
        prs$A_tau <- 0.01
        prs$B_tau <- 0.01
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
    JAGS_data <- c(JAGS_data, prs)
    ######################################################################################
    # write model
    model_name <- paste0("jags_mv_glmer", sample(1e06, 1), ".txt")
    cat(build_model(families, seq_along(families), colmns_HC, colmns_nHC, overdispersion,
                    JAGS_data$n_RE), file = file.path(con$working.directory, model_name))
    closeAllConnections()
    # parameters to save
    params <- paste0('betas', seq_len(n_outcomes))
    if (any(ind_gs <- sapply(families, function (x) x$family == "gaussian"))) {
        params <- c(params, paste0("sigma", which(ind_gs)))
    }
    params <- c(params, "inv.D", "b")
    inits <- function () {
        ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
                        rnorm, sd = 0.1)
        names(ints) <- paste0('betas', seq_len(n_outcomes))
        ints$u <- drop(matrix(rnorm(JAGS_data$n * JAGS_data$n_RE), JAGS_data$n,
                              JAGS_data$n_RE))
        ints$inv.D <- if (JAGS_data$n_RE > 1) diag(JAGS_data$n_RE) else 1
        if (any(ind_gs)) {
            nms <- which(ind_gs)
            taus <- rep(list(1), length(nms))
            names(taus) <- paste0("tau", nms)
            ints <- c(ints, taus)
        }
        ints
    }
    jags_fit <- jagsUI::jags(data = JAGS_data, inits = inits, parameters.to.save = params,
                     model.file = file.path(con$working.directory, model_name),
                     parallel = con$n.processors > 1, n.chains = con$n.chains,
                     n.adapt = con$n.adapt, n.iter = con$n.iter, n.burnin = con$n.burnin,
                     n.thin = con$n.thin, seed = con$seed, verbose = con$verbose)
    if (con$clear.model) {
        file.remove(file.path(con$working.directory, model_name))
    }
    out <- list(mcmc = jags_fit$sims.list, components = components, data = data,
                families = families, control = con, mcmc.info = jags_fit$mcmc.info,
                DIC = jags_fit$DIC, pD = jags_fit$pD, Rhat = jags_fit$Rhat,
                priors = prs)
    if (JAGS_data$n_RE == 1) {
        out$mcmc$inv.D <- array(out$mcmc$inv.D, c(length(out$mcmc$inv.D), 1, 1))
        out$mcmc$b <- array(out$mcmc$b, c(nrow(out$mcmc$b), ncol(out$mcmc$b), 1))
    }
    out$mcmc$D <- out$mcmc$inv.D
    for (i in seq_len(nrow(out$mcmc$betas1))) {
            out$mcmc$D[i, , ] <- solve(out$mcmc$D[i, , ])
    }
    # fix names
    pat <- paste0("^X", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
    Xnams <- lapply(components[grep(pat, names(components))], colnames)
    for (i in seq_along(Xnams)) {
        colnames(out$mcmc[[paste0("betas", i)]]) <- Xnams[[i]]
    }
    pat <- paste0("^Z", paste(rep("[0-9]", nchar(as.character(n_outcomes))), collapse = ""))
    Znams <- lapply(components[grep(pat, names(components))], colnames)
    Znams <- unlist(mapply(paste0, Znams, seq_len(n_outcomes), SIMPLIFY = FALSE))
    dimnames(out$mcmc$D) <- dimnames(out$mcmc$inv.D) <- list(NULL, Znams, Znams)
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
                DIC = object$DIC, pD = object$pD, call = object$call)
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
        print(x[[paste0("Outcome", i)]])
    }
    cat("\nMCMC summary:\n")
    tt <- x$mcmc.info$elapsed.mins
    cat("iterations:", x$control$n.iter, "\nadapt:", x$control$n.adapt,
        "\nburn-in:", x$control$n.burnin, "\nthinning:", x$control$n.thin,
        "\ntime:", if (tt > 60) round(tt/60, 1) else round(tt, 1),
        if (tt > 60) "hours" else "min")
    cat("\n")
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
