plot.rocJM <-
function (x, legend = FALSE, optimal.cutoff = c("", "qOverall", "Youden"),
                        ...) {
    plot(x$FP, x$TP, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity")
    abline(a = 0, b = 1, lty = 3)
    optimal.cutoff <- match.arg(optimal.cutoff)
    if (optimal.cutoff == "qOverall")
        abline(v = x$thrs[which.max(x$qOverall)], lty = 3, lwd = 2, col = 2)
    if (optimal.cutoff == "Youden")
        abline(v = x$thrs[which.max(x$TP - x$FP)], lty = 3, lwd = 2, col = 2)
    if (legend) {
        legend("bottomright", c(paste("At time:", round(x$Thoriz, 1), "\n"),
                                paste("\nUsing information up to time:", 
                                      round(x$Tstart, 1))), 
               bty = "n")
    }
    invisible()
}
