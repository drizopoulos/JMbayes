rmvt <-
function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
    } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
}
