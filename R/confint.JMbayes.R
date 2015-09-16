confint.JMbayes <-
function (object, parm = c("all", "Longitudinal", "Event"), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    parm <- match.arg(parm)
    cf <- switch(parm,
                 "Longitudinal" = fixef(object),
                 "Event" = fixef(object, "Event"),
                 "all" = {
                     cy <- fixef(object)
                     names(cy) <- paste("Y.", names(cy), sep = "")
                     ct <- fixef(object, "Event")
                     names(ct) <- paste("T.", names(ct), sep = "")
                     c(cy, ct)
                 })
    pnames <- names(cf)
    a <- (1 - 0.95)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(cf), 3L), 
                dimnames = list(names(cf), c(pct[1], "est.", pct[2])))
    civals <- switch(parm,
                     "Longitudinal" = object$CIs$betas,
                     "Event" = cbind(object$CIs$gammas, object$CIs$alphas, object$CIs$Dalphas, object$CIs$shapes),
                     "all" = cbind(object$CIs$betas, object$CIs$gammas, object$CIs$alphas, object$CIs$Dalphas, 
                                   object$CIs$shapes)
    )
    ci[, c(1,3)] <- t(civals)
    ci[, 2] <- cf
    ci
}
