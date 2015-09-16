update.JMbayes <-
function (object, ...) {
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component.\n")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        nams <- names(extras)
        existing <- !is.na(match(nams, names(call)))
        for (a in names(extras)[existing]) {
            call[[a]] <- extras[[a]]
        }
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
        if (all(nams %in% c("scales", "n.iter", "n.burnin", "n.adapt", "n.thin"))) {
            call <- c(as.list(call), list(init = extractInits(object)))
            call <- as.call(call)
        }
    } else {
        call <- c(as.list(call), list(init = extractInits(object)))
        call <- as.call(call)
    }
    eval(call, parent.frame())
}
