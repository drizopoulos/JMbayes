print.rocJM <-
function (x, digits = 4, ...) {
    cat("\n\tTime-dependent Sensitivity and Specificity for the Joint Model",  
            x$nameObject)
    cat("\n\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing information up to time: ", round(x$Tstart, digits), 
        " (", x$nr, " subjects still at risk)\n\n", sep = "")
    d <- data.frame("cut-off" = x$thrs, "SN" = x$TP, "SP" = 1 - x$FP, "qSN" = x$qSN, 
                    "qSP" = x$qSP, "qOverall" = x$qOverall, check.names = FALSE,
                    check.rows = FALSE)
    xx <- rep("", nrow(d))
    xx[ind <- which.max(d$qOverall)] <- "*"
    d[[" "]] <- xx
    nms <- c(as.character(seq(6, 96, 5)), row.names(d)[ind])
    d <- d[row.names(d) %in% nms, ]
    d <- d[!is.na(d$qOverall), ]
    row.names(d) <- 1:nrow(d)
    print(d)
    cat("\n")
    invisible(x)
}
