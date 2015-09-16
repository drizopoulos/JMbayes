trimAboveMean <-
function (x, trim = 0.1) {
    qx <- quantile(x, probs = 1 - trim, names = FALSE, na.rm = TRUE)
    x[x >= qx] <- as.numeric(NA)
    mean(x, na.rm = TRUE)
}
