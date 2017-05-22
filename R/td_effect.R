td <- function (x, df = NULL, knots = NULL, ord = 3) {
    if (is.null(df) && is.null(knots)) {
        stop("either 'df' or 'knots' need to be specified.\n")
    }
    if (is.null(knots) && !is.null(df)) {
        min_x <- min(x) - 0.1
        max_x <- max(x) + 0.1
        dx <- (max_x - min_x) / df
        knots <- seq(min_x - (ord-1) * dx, max_x + (ord-1) * dx, by = dx)
    }
    out <- splineDesign(knots, x, ord = ord)
    attr(out, 'knots') <- knots
    attr(out, "class") <- c("td", "basis", "matrix")
    out
}

makepredictcall.td <- function (var, call) {
    if (as.character(call)[1L] != "td") 
        return(call)
    at <- attributes(var)[c("knots", "ord")]
    xxx <- call[1L:2L]
    xxx[names(at)] <- at
    xxx
}