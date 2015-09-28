ModelMats <- function (time, ii, timeL = NULL) {
    GQsurv <- if (object$control$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(object$control$GQsurv.k)
    wk <- GQsurv$wk
    sk <- GQsurv$sk
    id.GK <- rep(ii, each = length(sk))
    if (!is.null(timeL)) {
        P <- c(time - timeL) / 2
        st <- P * sk +  c(time + timeL) / 2
    } else {
        P <- time / 2
        st <- P * (sk + 1)
    }
    #data.id2 <- data.id[id.GK, ]
    data.id2 <- newdata[id == ii, ]
    data.id2 <- data.id2[findInterval(st, obs.times[[ii]]), , drop = FALSE]
    data.id2[[timeVar]] <- pmax(st - lag, 0)
    kn <- object$control$knots
    W2s <- splineDesign(unlist(kn, use.names = FALSE), st, 
                        ord = object$control$ordSpline, outer.ok = TRUE)    
    out <- list(st = st, wk = rep(wk, length(P)), P = P, W2s = W2s, kn = kn, 
                idT = rep(seq_along(P), each = length(sk)))
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame.default(delete.response(TermsX), data = data.id2)
        mfZ <- model.frame.default(TermsZ, data = data.id2)
        out$Xs <- model.matrix.default(formYx, mfX)
        out$Zs <- model.matrix.default(formYz, mfZ)
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame.default(TermsX.extra, data = data.id2)
        mfZ.extra <- model.frame.default(TermsZ.extra, data = data.id2)
        out$Xs.extra <- model.matrix.default(extraForm$fixed, mfX.extra)
        out$Zs.extra <- model.matrix.default(extraForm$random, mfZ.extra)
    }
    if (estimateWeightFun) {
        P2 <- st / 2
        st2 <- outer(P2, sk + 1)
        id.GK2 <- rep(seq_len(nrow(data.id2)), each = length(sk))
        data.id3 <- data.id2[id.GK2, ]
        data.id3[[timeVar]] <- pmax(c(t(st2)) - lag, 0)
        mfX <- model.frame.default(delete.response(TermsX), data = data.id3)
        mfZ <- model.frame.default(TermsZ, data = data.id3)
        out$Xu <- model.matrix.default(formYx, mfX)
        out$Zu <- model.matrix.default(formYz, mfZ)
        out$P2 <- P2
        out$st2 <- st[id.GK2] - c(t(st2))
        out$wk2 <- rep(wk, length(P2))
        out$id.GK2 <- id.GK2
    }
    out
}
