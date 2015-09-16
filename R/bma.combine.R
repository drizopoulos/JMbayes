bma.combine <-
function (..., JMlis = NULL, weights = NULL) {
    lis <- if (is.null(JMlis) || !is.list(JMlis)) list(...) else JMlis 
    classes <- sapply(lis, class)
    if (is.matrix(classes))
        classes <- classes[2, ]
    if (!all(classes %in% c("survfit.JMbayes", "predict.JMbayes")))
        stop("bma.combine() only works for 'survfit.JMbayes' or 'predict.JMbayes' objects.\n")
    if (!all(classes == classes[1]))
        stop("all objects have to be of the same class.\n")
    K <- length(lis)
    if (all(classes == "survfit.JMbayes")) {    
        out <- lis[[1]]
        n <- length(out$summaries)
        if (is.null(weights))
            weights <- rep(1, K)
        for (i in 1:n) {
            out$summaries[[i]][, -1] <- out$summaries[[i]][, -1] * weights[1]
        }
        if (K > 1) {
            for (i in 1:n) {
                for (k in 2:K) {
                    out$summaries[[i]][, -1] <- out$summaries[[i]][, -1] + lis[[k]]$summaries[[i]][, -1] * weights[k]
                }
            }
        }
    } else {
        if (is.data.frame(lis[[1]])) {
            out <- lis[[1]]
            na.ind <- !is.na(out$se.fit)
            for (k in 2:K) {
                out$pred[na.ind] <- out$pred[na.ind] + lis[[k]][["pred"]][na.ind]
                out$low[na.ind] <- out$low[na.ind] + lis[[k]][["low"]][na.ind]
                out$upp[na.ind] <- out$upp[na.ind] + lis[[k]][["upp"]][na.ind]
            }
        } else {
            out <- lis[[1]]
            for (k in 2:K) {
                out$pred <- out$pred + lis[[k]][["pred"]]
                out$low <- out$low + lis[[k]][["low"]]
                out$upp <- out$upp + lis[[k]][["upp"]]
            }
        }
    }
    out
}
