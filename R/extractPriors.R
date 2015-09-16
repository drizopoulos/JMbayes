extractPriors <-
function (object, fitted = FALSE) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!fitted) {
        object$priors
    } else {
        pM <- object$postMeans
        list(priorMean.betas = pM$betas)
    }
}
