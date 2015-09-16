chol.transf <-
function (x) {
    if (any(is.na(x) | !is.finite(x))) 
        stop("NA or infinite values in 'x'.\n")
    if (is.matrix(x)) {
        k <- nrow(x)
        U <- chol(x)
        U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
        U[upper.tri(U, TRUE)]
    } else {
        nx <- length(x)
        k <- round((-1 + sqrt(1 + 8 * nx))/2)
        mat <- matrix(0, k, k)
        mat[upper.tri(mat, TRUE)] <- x
        mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
        res <- crossprod(mat)
        res
    }
}
