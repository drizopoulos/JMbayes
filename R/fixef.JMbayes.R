fixef.JMbayes <-
function (object, process = c("Longitudinal", "Event"), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        object$postMeans$betas
    } else {
        gammas <- object$postMeans$gammas
        out <- c(gammas, object$postMeans$alphas, object$postMeans$Dalphas, object$postMeans$shapes)
        if ((lag <- object$y$lag) > 0) {
            kk <- grep("Assoct", names(out), fixed = TRUE)
            names(out)[kk] <- paste(names(out)[kk], "(lag=", lag, ")", sep = "")
        }
        out
    }
}
