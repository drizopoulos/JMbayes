dwish <-
function (W, S, v, log = FALSE) {
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (!is.matrix(W)) 
        W <- matrix(W)
    k <- nrow(S)
    log.gammapart <- 0
    for (i in seq_len(k)) {
        log.gammapart <- log.gammapart + lgamma(0.5 * (v + 1 - i))
    }
    log.denom <- log.gammapart + (0.5 * v * k) * log(2) + (0.25 * k * (k - 1)) * log(pi)
    log.detS <- c(determinant(S, logarithm = TRUE)$modulus)
    log.detW <- c(determinant(W, logarithm = TRUE)$modulus)
    trace <- sum(diag(solve(S, W)))
    log.num <- - 0.5 * (v * log.detS - (v - k - 1) * log.detW + trace)
    if (log) 
        log.num - log.denom
    else 
        exp(log.num - log.denom)
}
