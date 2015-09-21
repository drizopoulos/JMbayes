dropAttr <- function (x, indBetas = NULL, iF = NULL) {
    if (is.null(x))
        return(NULL)
    if (is.matrix(x)) {
        d <- dim(x)
        x <- as.vector(x)
        dim(x) <- d
        if (!is.null(indBetas)) {
            if (is.null(iF)) {
                x <- x[, -indBetas, drop = FALSE]
            } else {
                x <- x[, !iF %in% indBetas, drop = FALSE]
            }
        }
        x
    } else if (is.list(x)) {
        lapply(x, unname)
    } else
        as.vector(x)
}
