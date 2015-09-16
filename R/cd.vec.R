cd.vec <-
function (x, f, ..., eps = .Machine$double.eps^0.25) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- eps * (abs(x) + eps)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + ex[i]
        x2[i] <- x[i] - ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f / diff.x
    }
    0.5 * (res + t(res))
}
