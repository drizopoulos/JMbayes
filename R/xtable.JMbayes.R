xtable.JMbayes <-
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, which = c("all", "Longitudinal", "Event"), varNames.Long = NULL, 
    varNames.Event = NULL, ...) {
    which <- match.arg(which)
    smobj <- summary(x)
    long <- smobj$'CoefTable-Long'
    long <- long[, - c(2:3)]
    if (!is.null(x$postMeans$sigma)) {
        sigma <- cbind("Value" = x$postMeans$sigma, "2.5%" = x$CIs$sigma[1], "97.5%" = x$CIs$sigma[2], "P" = x$Pvalues$sigma)
        rownames(sigma) <- "$\\sigma$"
    }
    D <- x$postMeans$D
    ind <- lower.tri(D, TRUE)
    D <- cbind("Value" = D[ind], "2.5%" = x$CIs$D[1, c(ind)], "97.5%" = x$CIs$D[2, c(ind)], "P" = x$Pvalues$D[c(ind)])
    Dat.long <- data.frame(rbind(long, sigma, D), check.names = FALSE)
    names(Dat.long) <- c("Value", "2.5\\%", "97.5\\%", "P")
    rnm.l <- row.names(Dat.long)
    if (!is.null(varNames.Long) && is.character(varNames.Long) && 
            length(varNames.Long) == length(rnm.l))
        row.names(Dat.long) <- varNames.Long
    Dat.long <- cbind(" " = row.names(Dat.long), Dat.long)
    event <- smobj$'CoefTable-Event'
    event <- event[, -c(2:3)]
    Dat.surv <- data.frame(event, check.names = FALSE)
    names(Dat.surv) <- c("Value", "2.5\\%", "97.5\\%", "P")
    rnm.e <- row.names(Dat.surv)
    if (!is.null(varNames.Event) && is.character(varNames.Event) && 
        length(varNames.Event) == length(rnm.e))
        row.names(Dat.surv) <- varNames.Event
    Dat.surv <- cbind(" " = row.names(Dat.surv), Dat.surv)
    nn <- nrow(Dat.long) - nrow(Dat.surv)
    Dat <- if (which == "all") {
        if (is.null(caption))
            caption <- paste("Parameter estimates and 95\\% credibility intervals",
                "under the joint modeling analysis.",
                "D[i, j] denote the $ij$-element of the covariance matrix",
                "for the random effects.")
        if (nn > 0) {
            dd <- Dat.surv[seq_len(nn), ]
            dd[] <- lapply(dd, function (x) {x[] <- NA; x})
            Dat.surv <- rbind(Dat.surv, dd)
        }
        if (nn < 0) {
            dd <- Dat.long[seq_len(abs(nn)), ]
            dd[] <- lapply(dd, function (x) {x[] <- NA; x})
            Dat.long <- rbind(Dat.long, dd)
        }
        align <- "llrrrrlrrrr"
        cbind(Dat.surv, Dat.long)
    } else if (which == "Longitudinal") {
        if (is.null(caption))
            caption <- paste("Parameter estimates and 95\\% credibility intervals",
                "under the joint modeling analysis for the longitudinal linear mixed effects submodel.",
                "D[i, j] denote the $ij$-element of the covariance matrix",
                "for the random effects.")
        align <- "llrrrr"
        Dat.long
    } else {
        if (is.null(caption))
            caption <- paste("Parameter estimates and 95\\% credibility intervals",
                "under the joint modeling analysis for the event time submodel.")
        align <- "llrrrr"
        Dat.surv
    }  
    if (!requireNamespace("xtable", quietly = TRUE)) 
        stop("'xtable' is required.\n")
    xtable::print.xtable(xtable::xtable(Dat, caption = caption, label = label, 
        align = align, digits = digits, display = display), 
        sanitize.text.function = function (x) x, include.rownames = FALSE, ...)
}
