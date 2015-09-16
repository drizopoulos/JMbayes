fd <-
function (x, f, ..., eps = .Machine$double.eps^0.25) {
    n <- length(x)
    res <- numeric(n)
    ex <- eps * (abs(x) + eps)
    f0 <- f(x, ...)
    for (i in seq_len(n)) {
        x1 <- x
        x1[i] <- x[i] + ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[i] <- diff.f / diff.x
    }
    res
}
