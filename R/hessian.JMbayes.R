hessian.JMbayes <-
function (object, thetas, numerDeriv = c("fd", "cd")) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (missing(thetas))
        thetas <- object$postMeans
    if (!is.list(thetas) || length(thetas) != length(object$postMeans[-3L]))
        stop("'thetas' must be a list with the model's parameters with the same structure as ",
            "'object$postMeans'.")    
    numerDeriv <- match.arg(numerDeriv)
    list.thetas <- thetas
    p <- ncol(list.thetas$D)
    list.thetas$D <- list.thetas$D[lower.tri(list.thetas$D, TRUE)]
    vec.thetas <- unlist(as.relistable(list.thetas))
    lL <- function (thetas) {
        tht <- relist(thetas, list.thetas)
        D <- matrix(0, p, p)
        D[lower.tri(D, TRUE)] <- tht$D
        tht$D <- D + t(D)
        diag(tht$D) <- 0.5 * diag(tht$D)
        - logLik(object, thetas = tht)
    }
    sc <- function (thetas) {
        if (numerDeriv == "fd") fd(thetas, lL) else cd(thetas, lL)
    }
    if (numerDeriv == "fd") fd.vec(vec.thetas, sc) else cd.vec(vec.thetas, sc)
}
