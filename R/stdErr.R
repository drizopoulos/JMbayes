stdErr <-
function (x) {
    x <- as.matrix(x)
    vars <- apply(x, 2L, var)
    ess <- effectiveSize(x)
    sqrt(vars / ess)
}
