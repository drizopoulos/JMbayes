computeP <-
function (x) {
    above <- mean(x >= 0)
    below <- mean(x < 0)
    2 * min(above, below)
}
