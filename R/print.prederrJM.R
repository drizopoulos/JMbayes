print.prederrJM <-
function (x, digits = 4, ...) {
    if (!inherits(x, "prederrJM"))
        stop("Use only with 'prederrJM' objects.\n")
    if (x$class == "JMbayes" || x$class == "jointModel" || x$class == "mvJMbayes")
        cat("\nPrediction Error for the Joint Model", x$nameObject)
    else
        cat("\nPrediction Error for the Cox model", x$nameObject)
    cat("\n\nEstimated prediction error:", round(x$prederr, digits))
    if (!x$interval) {
        cat("\nAt time:", round(x$Thoriz, digits))
    } else {
        cat("\nIn the time interval: [", round(x$Tstart, digits), 
            ", ", round(x$Thoriz, digits), "]", sep = "")
    }
    cat("\nUsing information up to time: ", round(x$Tstart, digits), " (", x$nr, " subjects still at risk)", sep = "")
    cat("\nLoss function:", if (is.function(x$lossFun)) "user-defined function" else x$lossFun)
    cat("\n\n")
    invisible(x)
}
