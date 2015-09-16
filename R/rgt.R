rgt <-
function (n, mu = 0, sigma = 1, df = stop("no df arg")) {
    mu + sigma * rt(n, df = df)
}
