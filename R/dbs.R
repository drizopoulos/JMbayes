dbs <-
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), eps = 1e-03) {
    bs.x <- if (is.null(knots)) {
        splines::bs(x, df = df, intercept = intercept, Boundary.knots = Boundary.knots)
    } else {
        splines::bs(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
    } 
    kn <- attr(bs.x, "knots")
    Bkn <- attr(bs.x, "Boundary.knots")
    ex <- pmax(abs(x), 1)
    x1 <- x + eps * ex
    x2 <- x - eps * ex
    bs.xeps1 <- suppressWarnings(splines::bs(x1, knots = kn, Boundary.knots = Bkn, intercept = intercept))
    bs.xeps2 <- suppressWarnings(splines::bs(x2, knots = kn, Boundary.knots = Bkn, intercept = intercept))
    out <- (bs.xeps1 - bs.xeps2) / c(x1 - x2)
    attr(out, "eps") <- eps
    attr(out, "class") <- c("dbs", "basis", "matrix")
    out
}
