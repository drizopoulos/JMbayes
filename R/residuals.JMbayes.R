residuals.JMbayes <-
function (object, process = c("Longitudinal", "longitudinal", "Event", "event"), 
                               type = c("Marginal", "marginal", "Subject", "subject", 
                                        "Martingale", "martingale", "nullMartingale", "nullmartingale"),
                               standardized = FALSE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal" || process == "longitudinal") {
        if (type %in% c("Martingale", "martingale", "nullMartingale", "nullmartingale"))
            stop("invalid 'type'; it can only be 'marginal' or 'subject'.\n")
        y <- object$y$y
        fits <- fitted(object, process = "Longitudinal", type = type)
        res <- y - fits
        if (standardized && type %in% c("Subject", "subject"))
            res <- res / object$postMeans$sigma
        res
    } else {
        ni <- tapply(object$y$id, object$y$id, length)
        events <- rep(object$y$event, ni)
        events <- ave(events, object$y$id, FUN = function (x) c(rep(0, length(x) - 1L), x[1L]))                                            
        fits <- fitted(object, process = "Event", 
                       nullY = type %in% c("nullMartingale", "nullmartingale"))
        events - fits
    }
}
