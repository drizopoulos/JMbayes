print.summary.JMbayes <-
function (x, digits = max(4, getOption("digits") - 4),
                                   printKnots = FALSE, ...) {
    if (!inherits(x, "summary.JMbayes"))
        stop("Use only with 'summary.JMbayes' objects.\n")
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("Data Descriptives:\n")
    pcEv <- round(100 * sum(x$d) / x$n, 1)
    cat("Longitudinal Process\t\tEvent Process")
    cat("\nNumber of Observations: ", x$N, "\tNumber of Events: ",
        sum(x$d), " (", pcEv, "%)", sep = "")
    cat("\nNumber of subjects:", x$n)
    cat("\n\nJoint Model Summary:")
    if (x$densLongCheck) {
        if (is.null(x$df.RE)) {
            cat("\nLongitudinal Process: Linear mixed-effects model")
        } else {
            cat("\nLongitudinal Process: Linear mixed-effects model with Student's-t(df=", 
                x$df.RE, ") random effects", sep = "")        
        }
    } else {
        if (is.null(x$df.RE)) {
            cat("\nLongitudinal Process: user-defined mixed model")
        } else {
            cat("\nLongitudinal Process: user-defined mixed model with Student's-t(df=", 
                x$df.RE, ") random effects", sep = "")        
        }
    }
    cat("\nEvent Process: ")
    xx <- if (length(x$control$knots) == 1L) {
        kk <- round(unique(x$control$knots[[1]]), 1)
        paste(kk[-c(1, length(kk))], collapse = ", ")
    } else {
        paste(names(x$control$knots), sapply(x$control$knots, function (k) {
            kk <- round(unique(k), 1)
            paste(kk[-c(1, length(kk))], collapse = ", ")
        }), sep = ": ", collapse = "\n\t\t")
    }
    ttE <- if (x$baseHaz == "P-splines") "penalized-spline-approximated" else "spline-approximated"
    if (printKnots)
        cat("Relative risk model with ", ttE, " baseline risk function (knots at: ", 
            xx, ")\n", sep = "")
    else
        cat("Relative risk model with", ttE, "\n\t\tbaseline risk function\n")
    if (x$estimateWeightFun) {
        cat("Parameterization: weighted cumulative effect\n\n")
    } else {
        cat("Parameterization:", switch(x$param, "td-value" = "Time-dependent value",
                                        "td-extra" = "Time-dependent extra term", 
                                        "td-both" = "Time-dependent value + time-dependent extra term",
                                        "shared-betasRE" = "shared subject-specific coefficients",
                                        "shared-RE" = "shared random effects"), "\n\n")
    }
    if (!is.null(x$DIC)){
        model.sum <- data.frame(LPML = x$LPML, DIC = x$DIC, pD = x$pD, row.names = "")
        print(model.sum)
    }
    cat("\nVariance Components:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- as.data.frame(round(rbind(sds, "Residual" = x$sigma), digits))
            names(dat) <- "StdDev"
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- rbind(mat, c(x$sigma, rep(0, ncz - 1)))
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- if (!is.null(x$sigma)) 
                c(dimnames(D)[[1]], "Residual") else c(dimnames(D)[[1]], " ")
        }
    } else {
        dat <- data.frame("StdDev" = c(sds, x$sigma), 
                          row.names = if (!is.null(x$sigma)) c(rownames(D), "Residual") else rownames(D),
                          check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    cat("\nCoefficients:")
    cat("\nLongitudinal Process\n")
    out <- as.data.frame(round(x$"CoefTable-Long", digits))
    out$P <- format.pval2(out$P, digits = digits, eps = 1e-03)
    print(out)
    cat("\nEvent Process\n")
    out <- as.data.frame(round(x$"CoefTable-Event", digits))
    out$P <- format.pval2(out$P, digits = digits, eps = 1e-03)
    print(out)
    cat("\nMCMC summary:\n")
    tt <- x$time["elapsed"]/60
    cat("iterations:", x$control$n.iter, "\nadapt:", x$control$n.adapt, 
        "\nburn-in:", x$control$n.burnin, "\nthinning:", x$control$n.thin,
        "\ntime:", if (tt > 60) round(tt/60, 1) else round(tt, 1), 
        if (tt > 60) "hours" else "min")
    cat("\n")
    invisible(x)
}
