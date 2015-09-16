print.dynCJM <-
function (x, digits = 4, ...) {
    if (!inherits(x, "dynCJM"))
        stop("Use only with 'dynCJM' objects.\n")
    if (x$classObject == "JMbayes")
        cat("\n\tDynamic Discrimination Index for the Joint Model",  x$nameObject)
    else
        cat("\n\tDynamic Discrimination Index for the Cox Model",  x$nameObject)
    cat("\n\nEstimated dynC:", round(x$dynC, digits))
    cat("\nIn the time interval: [0, ", round(x$t.max, digits), "]", sep = "")
    cat("\nLength of time interval:", round(x$Dt, digits))
    cat("\n\n")
    invisible(x)
}
