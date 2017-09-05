print.survfit.JMbayes <-
function (x, ...) {
    if (!is.null(x$success.rate)) {
        cat("\nPrediction of Conditional Probabilities of Event\n\tbased on", 
            nrow(x$success.rate), "Monte Carlo samples\n\n")
    } else {
        cat("\nPrediction of Conditional Probabilities for Events\n\n")
    }
    f <- function (d, t) {
        dd <- d[1L, , drop = FALSE]
        dd[1L, ] <- c(as.vector(t), rep(1, ncol(dd) - 1))
        round(rbind(dd, d), 4L)
    }
    print(mapply(f, x$summaries, x$last.time, SIMPLIFY = FALSE))
    invisible(x)
}
