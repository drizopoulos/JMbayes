plot.JMbayes <-
function (x, which = c("trace", "autocorr", "density", "CPO", "weightFun"), 
                          param = c("betas", "sigma", "D", "gammas", "alphas", "Dalphas", 
                                    "shapes", "Bs.gammas", "tauBs"), 
                          ask = TRUE, max.t = NULL, from = 0, ...) {
    if (!inherits(x, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    which <- match.arg(which)
    if (which %in% c("trace", "density", "autocorr")) {
        param <- match.arg(param, several.ok = TRUE)
        if (any(param == "D")) {
            keepD <- c(lower.tri(x$postMeans$D, TRUE))
            x$mcmc$D <- x$mcmc$D[, keepD]
        }
        pp <- do.call(cbind, x$mcmc[param])
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
    } else if (which == "CPO") {
        n <- length(x$CPO)
        matplot(matrix(1:n, 2, n, TRUE), rbind(rep(0, n), x$CPO), type = "l",
                xlab = "Subjects", ylab = "CPO", main = deparse(substitute(x)), ...)
    } else {
        if (!x$estimateWeightFun)
            stop("\nweight function has not been estimated.")
        weightFun <- x$Funs$weightFun
        if (is.null(max.t))
            max.t <- quantile(x$y$Time, 0.25)
        xx <- seq(1e-03, max.t, length.out = 101)
        yy <- weightFun(xx, x$postMeans$shapes, max(x$y$Time))
        plot(xx, yy, type = "l", xlab = "t - s", ylab = "weights", 
             main = "Estimated Weight Function")
    }
    invisible()
}
