check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z == z[1])))
}
