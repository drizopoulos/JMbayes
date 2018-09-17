summary.mvJMbayes <- function (object, weighted = FALSE, include.baselineHazard = FALSE, ...) {
    families <- object$model_info$families
    n_outcomes <- length(families)
    components <- object$model_info$mvglmer_components
    extract_components <- function (nam) {
        components[grep(nam, names(components), fixed = TRUE)]
    }
    respVars <- unlist(extract_components("respVar"), use.names = FALSE)
    descrpt <- data.frame(" " = unlist(extract_components("N"), use.names = FALSE),
                          row.names = respVars, check.rows = FALSE, check.names = FALSE)
    postMeans <- if (weighted && !is.null(object$statistics$postwMeans)) {
        object$statistics$postwMeans 
    } else {
        object$statistics$postMeans
    }
    CIs <- if (weighted && !is.null(object$statistics$wCIs)) {
        object$statistics$wCIs
    } else {
        object$statistics$CIs
    }
    StDev <- if (weighted && !is.null(object$statistics$wStDev)) {
        object$statistics$wStDev
    } else {
        object$statistics$StDev
    }
    out <- list(n = components$n1, descrpt = descrpt, D = postMeans$D,
                families = families, respVars = respVars, 
                events = object$model_info$coxph_components$event,
                control = object$control, mcmc_info = object$mcmc_info, call = object$call)
    tab_f <- function (name, is_sigma = FALSE) {
        is_mat <- is.matrix(object$statistics$CIs[[name]])
        data.frame("PostMean" = postMeans[[name]],
                   "StDev" = StDev[[name]],
                   "StErr"= object$statistics$StErr[[name]],
                   "2.5%" = if (is_mat) CIs[[name]][1, ] else CIs[[name]][1],
                   "97.5%" = if (is_mat) CIs[[name]][2, ] else CIs[[name]][2],
                   "P" = object$statistics$Pvalues[[name]],
                   row.names = if (is_sigma) "sigma" 
                   else names(object$statistics$postMeans[[name]]),
                   check.names = FALSE)
    }
    for (i in seq_len(n_outcomes)) {
        out[[paste0("Outcome", i)]] <- tab_f(paste0("betas", i))
        if (families[[i]][["family"]] == "gaussian") {
            D <- tab_f(paste0("sigma", i), TRUE)
            out[[paste0("Outcome", i)]] <- rbind(out[[paste0("Outcome", i)]], D)
        }
    }
    out$Survival <- do.call(rbind, list(tab_f("gammas"), tab_f("alphas"),
                                        if (include.baselineHazard) tab_f("Bs_gammas")))
    class(out) <- "summary.mvJMbayes"
    out
}

print.summary.mvJMbayes <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("Data Descriptives:")
    cat("\nNumber of Groups: ", x$n, "\t\tNumber of events: ", sum(x$event == 1),
        " (", round(100 * mean(x$event == 1), 1), "%)", sep = "")
    cat("\nNumber of Observations:")
    obs <- x$descrpt
    for (i in 1:nrow(obs)) {
        cat("\n  ", row.names(obs)[i], ": ", obs[[1]][i], sep = "")
    }
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
            row.names(dat) <- abbreviate(c(dimnames(D)[[1]]))
        }
    } else {
        dat <- data.frame("StdDev" = c(sds, x$sigma),
                          row.names = if (!is.null(x$sigma)) c(rownames(D), "Residual") else rownames(D),
                          check.rows = FALSE, check.names = FALSE)
    }
    print(dat, digits = digits)
    cat("\nSurvival Outcome:\n")
    print(round(x[["Survival"]], digits))
    n_outcomes <- length(x$families)
    for (i in seq_len(n_outcomes)) {
        cat("\nLongitudinal Outcome: ", x$respVars[i], " (family = ",
            x$families[[i]][["family"]], ", link = ", x$families[[i]][["link"]],
            ")", "\n", sep = "")
        xx <- round(x[[paste0("Outcome", i)]], digits)
        rnams <- row.names(xx)
        if (any(offend <- nchar(rnams) > 20))
            row.names(xx)[offend] <- abbreviate(rnams[offend])
        print(xx)
    }
    cat("\nMCMC summary:\n")
    tt <- x$mcmc_info$elapsed_mins
    cat("iterations:", x$mcmc_info$n_iter,
        "\nburn-in:", x$mcmc_info$n_burnin, "\nthinning:", x$mcmc_info$n_thin,
        "\ntime:", if (tt > 60) round(tt/60, 1) else round(tt, 1),
        if (tt > 60) "hours" else "min")
    cat("\n")
    invisible(x)
}

plot.mvJMbayes <- function (x, which = c("trace", "autocorr", "density", "tv_effect"),
                            param = c("betas", "sigma", "D", "gammas",
                                      "alphas", "Bs_gammas", "tau_Bs_gammas"),
                            ask = TRUE, max.t = NULL, from = 0, shade_CI = TRUE, 
                            col = 1, col_CI = "lightgrey", return_values = FALSE, ...) {
    if (!inherits(x, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    which <- match.arg(which)
    if (which %in% c("trace", "density", "autocorr")) {
        param <- match.arg(param, several.ok = TRUE)
        if (any(param == "D")) {
            keepD <- lower.tri(x$statistics$postMeans$D, TRUE)
            x$mcmc$D <- t(apply(x$mcmc$D, 1, c))[, c(keepD)]
            dnams <- which(keepD, arr.ind = TRUE)
            colnames(x$mcmc$D) <- paste0("D[", dnams[, 1], ", ", dnams[, 2], "]")
        }
        if (any(param == "tau_Bs_gammas")) {
            if (!x$model_info$multiState) {
                colnames(x$mcmc$tau_Bs_gammas) <- "tau_Bs_gammas"
            } else {
                colnames(x$mcmc$tau_Bs_gammas) <- c(paste0("tau_Bs_Gammas", "(trans = ", seq_len(ncol(x$mcmc$tau_Bs_gammas)), ")"))
            }
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
        values <- NULL
    }
    if (which == "tv_effect") {
        Data_surv <- x$model_info$coxph_components$data
        attr(Data_surv, "terms") <- NULL
        Time <- Data_surv[[1]][, 1]
        TimeVar <- all.vars(x$model_info$coxph_components$Terms)[1L]
        Data_surv[[TimeVar]] <- x$model_info$coxph_components$Time
        xx <- seq(min(Time), max(Time), length.out = 100)
        td_cols <- x$mcmc_info$priors$td_cols
        extract_alphas <- function (td_cols, alphas) alphas[td_cols]
        td_alphas_est <- lapply(td_cols, extract_alphas, alphas = x$statistics$postMeans$alphas)
        td_alphas_low <- lapply(td_cols, extract_alphas, alphas = x$statistics$CIs$alphas[1, ])
        td_alphas_upp <- lapply(td_cols, extract_alphas, alphas = x$statistics$CIs$alphas[2, ])
        exract_terms <- function (x, data = data) {
            environment(x) <- parent.frame()#.GlobalEnv
            terms.formula(x, data = data)
        }
        TD_terms <- lapply(x$model_info$Interactions, exract_terms, data = Data_surv)
        TD_frames <- lapply(TD_terms, model.frame.default, data = Data_surv)
        TD_terms <- lapply(TD_frames, terms)
        pred_data <- data.frame(V1 = xx)
        pred_data[[TimeVar]] <- pred_data$V1
        TD_mats <- lapply(TD_terms, model.matrix.default, data = pred_data)
        plot_fun <- function (est, low, upp, M, xx) {
            Fits <- cbind(est = M %*% est, low = M %*% low, upp = M %*% upp)
            if (return_values) {
                out <- cbind(times = xx, Fits)
                colnames(out) <- c("times", "est", "low", "upp")
                out
            } else {
                matplot(xx, Fits, type = "n", xlab = "Time", ylab = "Time-varying Effect", 
                        col = col)
                if (shade_CI) {
                    polygon(c(xx, rev(xx)), c(Fits[, 2], rev(Fits[, 3])), col = col_CI,
                            border = NA)
                }
                matlines(xx, Fits, col = 1, lty = c(1, 2, 2))
            }
        }
        values <- mapply(plot_fun, td_alphas_est, td_alphas_low, td_alphas_upp, TD_mats, 
                         MoreArgs = list(xx = xx), SIMPLIFY = FALSE)
    }
    invisible(values)
}

update.mvJMbayes <- function (object, ...) {
    call <- object$call
    if (is.null(call))
        stop("need an object with call component.\n")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        nams <- names(extras)
        existing <- !is.na(match(nams, names(call)))
        for (a in names(extras)[existing]) {
            call[[a]] <- extras[[a]]
        }
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    } else {
        call <- c(as.list(call), list(init = extractInits(object)))
        call <- as.call(call)
    }
    eval(call, parent.frame())
}

fixef.mvJMbayes <- function (object, process = c("Longitudinal", "Event"), ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        comps <- object$model_info$mvglmer_components
        nams_outcomes <- unlist(comps[grep("respVar", names(comps), fixed = TRUE)], 
                                use.names = FALSE)
        pMeans <- object$statistics$postMeans
        betas <- pMeans[grep("betas", names(pMeans), fixed = TRUE)]
        names(betas) <- nams_outcomes
        betas
    } else {
        c(object$statistics$postMeans$gammas, 
          object$statistics$postMeans$alphas)
    }
}

ranef.mvJMbayes <- function (object, as_list = FALSE, ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    out <- object$statistics$postMeans$b
    dimnames(out) <- list(names(object$model_info$coxph_components$Time), 
                          colnames(object$statistics$postMeans$D))
    if (as_list) {
        components <- object$model_info$mvglmer_components
        ncols <- components[grep("ncz", names(components), fixed = TRUE)]
        RE_inds <- object$model_info$RE_inds
        out <- lapply(RE_inds, function (i) out[, i, drop = FALSE])
    }
    out
}

fitted.mvJMbayes <- function (object, process = c("Longitudinal", "longitudinal", "Event", 
                                                  "event"), 
                              type = c("Marginal", "marginal", "Subject", "subject"), ...) {
    if (!inherits(object, "mvJMbayes")) 
        stop("Use only with 'mvJMbayes' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal" || process == "longitudinal") {
        components <- object$model_info$mvglmer_components
        families <- object$model_info$families
        n_outcomes <- length(families)
        seq_n_outcomes <- seq_len(n_outcomes)
        X <- components[paste0("X", seq_n_outcomes)]
        fitY <- mapply("%*%", X, fixef(object), SIMPLIFY = FALSE)
        names(fitY) <- row.names(object$Data$data)
        if (type == "Subject" || type == "subject") {
            ids <- components[paste0("id", seq_n_outcomes)]
            Zs <- components[paste0("Z", seq_n_outcomes)]
            bs <- ranef(object, as_list = TRUE)
            Zb_fun <- function (Z, b, id) rowSums(Z * b[id, , drop = FALSE])
            Zbs <- mapply(Zb_fun, Zs, bs, ids, SIMPLIFY = FALSE)
            fitY <- mapply("+", fitY, Zbs, SIMPLIFY = FALSE)
        }
        fitY
    }
}