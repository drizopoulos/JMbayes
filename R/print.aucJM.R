print.aucJM <-
function (x, digits = 4, ...) {
    if (!inherits(x, "aucJM"))
        stop("Use only with 'aucJM' objects.\n")
    if (x$class == "JMbayes") 
        cat("\n\tTime-dependent AUC for the Joint Model",  x$nameObject)
    else
        cat("\n\tTime-dependent AUC for the Cox Model",  x$nameObject)
    cat("\n\nEstimated AUC:", round(x$auc, digits))
        cat("\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing information up to time: ", round(x$Tstart, digits), 
        " (", x$nr, " subjects still at risk)", sep = "")
    cat("\n\n")
    invisible(x)
}
