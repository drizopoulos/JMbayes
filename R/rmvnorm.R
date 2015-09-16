rmvnorm <-
function (n, mu = NULL, Sigma) {
    if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
    } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
    }
    p <- length(ev)
    X <- tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), matrix(rnorm(n * p), n))
    if (!is.null(mu)) 
        X <- drop(mu) + X
    X <- if (n == 1L) drop(X) else t.default(X)
    X
}
