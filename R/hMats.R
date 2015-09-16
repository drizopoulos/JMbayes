hMats <-
function (Time) {
    W2 <- splineDesign(unlist(object$control$knots, use.names = FALSE), Time, 
                       ord = object$control$ordSpline, outer.ok = TRUE)
    data.id2 <- data.id
    data.id2[[timeVar]] <- pmax(Time - lag, 0)
    out <- list(W2 = W2, data = data.id2)
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame.default(delete.response(TermsX), data = data.id)
        mfZ <- model.frame.default(TermsZ, data = data.id)
        out$Xtime <- model.matrix.default(formYx, mfX)
        out$Ztime <- model.matrix.default(formYz, mfZ)
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame.default(TermsX.extra, data = data.id)
        mfZ.extra <- model.frame.default(TermsZ.extra, data = data.id)
        out$Xtime.extra <- model.matrix.default(extraForm$fixed, mfX.extra)
        out$Ztime.extra <- model.matrix.default(extraForm$random, mfZ.extra)
    }
    if (estimateWeightFun) {
        GQsurv <- if (object$control$GQsurv == "GaussKronrod") gaussKronrod() 
            else gaussLegendre(object$control$GQsurv.k)
        wk <- GQsurv$wk
        sk <- GQsurv$sk
        P <- Time / 2
        st <- outer(P, sk + 1)
        id.GK <- rep(seq_len(nrow(data.id2)), each = length(sk))
        data.id3 <- data.id2[id.GK, ]
        data.id3[[timeVar]] <- pmax(c(t(st)) - lag, 0)
        mfX <- model.frame.default(delete.response(TermsX), data = data.id3)
        mfZ <- model.frame.default(TermsZ, data = data.id3)
        out$Xu <- model.matrix.default(formYx, mfX)
        out$Zu <- model.matrix.default(formYz, mfZ)
        out$P <- P
        out$st <- Time[id.GK] - c(t(st))
        out$wk <- rep(wk, length(P))
        out$id.GK <- id.GK
        out$data <- data.id3
    }
    out
}
