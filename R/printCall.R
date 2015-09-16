printCall <-
function (call) {
    d <- deparse(call)
    if (length(d) <= 3) {
        paste(d, sep = "\n", collapse = "\n")
    } else {
        d <- d[1:3]
        d[3] <- paste0(d[3], "...")
        paste(d, sep = "\n", collapse = "\n")
    }
}
