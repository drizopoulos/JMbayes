rowMedians <-
function (mat) {
    apply(mat, 1, median, na.rm = TRUE)
}
