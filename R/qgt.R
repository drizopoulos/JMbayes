qgt <-
function (p, mu = 0, sigma = 1, df = stop("no df arg")) {
    mu + sigma * qt(p, df = df)
}
