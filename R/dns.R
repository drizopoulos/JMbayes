dns <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), eps = 1e-03) {
    ns.x <- if (is.null(knots)) {
        ns(x, df = df, intercept = intercept, Boundary.knots = Boundary.knots)
    } else {
        ns(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
    } 
    kn <- attr(ns.x, "knots")
    Bkn <- attr(ns.x, "Boundary.knots")
    ex <- pmax(abs(x), 1)
    x1 <- x + eps * ex
    x2 <- x - eps * ex
    ns.xeps1 <- ns(x1, knots = kn, Boundary.knots = Bkn, intercept = intercept)
    ns.xeps2 <- ns(x2, knots = kn, Boundary.knots = Bkn, intercept = intercept)
    out <- (ns.xeps1 - ns.xeps2) / c(x1 - x2)
    attr(out, "eps") <- eps
    attr(out, "class") <- c("dns", "basis", "matrix")
    out
}
