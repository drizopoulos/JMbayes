print.JMbayes <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
    if (!inherits(x, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("Variance Components:\n")
    D <- x$postMeans$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- round(c(D), digits)
            names(dat) <- rownames(D)
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(mat)[-1], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- dimnames(D)[[1]]
        }
    } else {
        dat <- data.frame("StdDev" = sds, row.names = rownames(D), 
                          check.rows = FALSE, check.names = FALSE)
    }
    lis <- list("Random-Effects" = dat, "Residual Std. Err." = x$postMeans$sigma)
    print(lapply(lis, function (x) if (!is.numeric(x)) x else round(x, digits = digits)))
    cat("\nCoefficients:\n")
    gammas <- c(x$postMeans$gammas, x$postMeans$alphas, x$postMeans$Dalphas, x$postMeans$shapes, x$postMeans$Bs.gammas)
    if ((lag <- x$y$lag) > 0) {
        kk <- grep("Assoct", names(gammas), fixed = TRUE)
        names(gammas)[kk] <- paste(names(gammas)[kk], "(lag=", lag, ")", sep = "")
    }
    print(lapply(list("Longitudinal Process" = x$postMeans$betas, "Event Process" = gammas), 
                 round, digits = digits))
    if (!is.null(x$DIC))
        cat("\nDIC:", x$DIC)
    cat("\n\n")
    invisible(x)
}
