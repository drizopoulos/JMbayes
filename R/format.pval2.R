format.pval2 <-
function (pv, digits = max(1L, getOption("digits") - 2L), eps = .Machine$double.eps, 
          na.form = "NA", ...) {
    if ((has.na <- any(ina <- is.na(pv)))) 
        pv <- pv[!ina]
    r <- character(length(is0 <- pv < eps))
    if (any(!is0)) {
        rr <- pv <- pv[!is0]
        expo <- floor(log10(ifelse(pv > 0, pv, 1e-50)))
        fixp <- expo >= -3 | (expo == -4 & digits > 1)
        if (any(fixp)) 
            rr[fixp] <- format(pv[fixp], digits = digits, ...)
        if (any(!fixp)) 
            rr[!fixp] <- format(pv[!fixp], digits = digits, ...)
        r[!is0] <- rr
    }
    if (any(is0)) {
        digits <- max(1L, digits - 2L)
        if (any(!is0)) {
            nc <- max(nchar(rr, type = "w"))
            if (digits > 1L && digits + 6L > nc) 
                digits <- max(1L, nc - 7L)
            sep <- if (digits == 1L && nc <= 6L) 
                ""
            else " "
        } else sep <- if (digits == 1) 
            ""
        else " "
        r[is0] <- paste("<", format(eps, digits = digits, ...), 
                        sep = "")
    }
    if (has.na) {
        rok <- r
        r <- character(length(ina))
        r[!ina] <- rok
        r[ina] <- na.form
    }
    r
}
