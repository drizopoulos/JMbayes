dgt <-
function (x, mu = 0, sigma = 1, df = stop("no df arg"), log = FALSE) {
    if (log)
        dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)
    else
         dt((x - mu) / sigma, df = df) / sigma
}
