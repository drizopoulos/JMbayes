dmvnorm <-
function (x, mu, Sigma = NULL, invSigma = NULL, log = FALSE, prop = TRUE) {
    if (!is.matrix(x))
        x <- rbind(x)
    if (is.matrix(mu)) {
        if (nrow(mu) != nrow(x))
            stop("incorrect dimensions for 'mu'.")
        p <- ncol(mu)
    } else {
        p <- length(mu)
        mu <- rep(mu, each = nrow(x))
    }
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
        invSigma <- evec %*% (t(evec) / ev)
        if (!prop)
            logdetSigma <- sum(log(ev))
    } else {
        if (!prop)
            logdetSigma <- - determinant(as.matrix(invSigma))$modulus
    }
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
    if (!prop)
        fact <- - 0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) {
        if (!prop) as.vector(fact - quad) else as.vector(- quad)
    } else {
        if (!prop) as.vector(exp(fact - quad)) else as.vector(exp(- quad))
    }
}
