pgt <-
function (q, mu = 0, sigma = 1, df = stop("no df arg")) {
    pt((q - mu) / sigma, df = df)
}
