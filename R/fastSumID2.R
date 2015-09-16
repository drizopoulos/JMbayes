fastSumID2 <-
function (x, group) {
    y <- c(0, cumsum(x)[group])
    out <- y[-1L] - y[-length(y)]
    if (anyNA(out))
        out[is.na(out)] <- -Inf
    out
}
