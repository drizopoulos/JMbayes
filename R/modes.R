modes <-
function (y) {
    test <- try(d <- density(y, bw = "nrd", adjust = 3, n = 1000), silent = TRUE)
    if (!inherits(test, "try-error")) d$x[which.max(d$y)] else NA
}
