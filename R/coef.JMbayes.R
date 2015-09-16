coef.JMbayes <-
function (object, process = c("Longitudinal", "Event"), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        betas <- object$postMeans$betas
        out <- matrix(betas, nrow = length(object$y$Time), ncol = length(betas), byrow = TRUE)
        colnames(out) <- names(betas)
        EB <- ranef(object)
        out[, colnames(EB)] <- out[, colnames(EB)] + EB
        out
    } else {
        out <- c(object$postMeans$gammas, object$postMeans$alphas, object$postMeans$Dalphas, object$postMeans$shapes)
          if ((lag <- object$y$lag) > 0) {
            kk <- grep("Assoct", names(out), fixed = TRUE)
            names(out)[kk] <- paste(names(out)[kk], "(lag=", lag, ")", sep = "")
        }
        out
    }
}
