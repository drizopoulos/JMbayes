dgbeta <-
function (x, shape1, shape2, d = 1, log = FALSE) {
    den <- lbeta(shape1, shape2) + (shape1 + shape2 - 1) * log(d)
    num <- (shape1 - 1) * log(x) + (shape2 - 1) * log(d - x)
    if (log) num - den else exp(num - den)
}
