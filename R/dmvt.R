dmvt <-
function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
        stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
        x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
        stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
        if (is.list(Sigma)) {
            ev <- Sigma$values
            evec <- Sigma$vectors
        } else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
        }
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            stop("'Sigma' is not positive definite")
        invSigma <- evec %*% (t(evec)/ev)
        if (!prop)
            logdetSigma <- sum(log(ev))
    } else {
        if (!prop)
            logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
        fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
                    log(df)) + logdetSigma)
    if (log) {
        if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
        if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
}
