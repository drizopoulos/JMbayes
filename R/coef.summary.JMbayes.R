coef.summary.JMbayes <-
function (object, ...) {
    if (!inherits(object, "summary.JMbayes"))
        stop("Use only with 'summary.JMbayes' objects.\n")
    coefsY <- object$'CoefTable-Long'
    coefsT <- object$'CoefTable-Event'
    list("Longitudinal" = coefsY, "Event" = coefsT)
}
