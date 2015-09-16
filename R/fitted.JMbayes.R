fitted.JMbayes <-
function (object, process = c("Longitudinal", "longitudinal", "Event", "event"), 
                            type = c("Marginal", "marginal", "Subject", "subject"), nullY = FALSE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal" || process == "longitudinal") {
        fitY <- c(object$x$X %*% object$postMeans$betas)
        names(fitY) <- row.names(object$Data$data)
        if (type == "Subject" || type == "subject")
            fitY <- fitY + rowSums(object$x$Z * ranef(object)[object$y$id, ])
        fitY
    } else {
        Data <- object$Data
        Funs <- object$Funs
        Forms <- object$Forms
        timeVar <- object$timeVar
        param <- object$param
        indFixed <- Forms$extraForm$indFixed
        indRandom <- Forms$extraForm$indRandom
        lag <- object$y$lag
        TermsX <- object$Terms$termsYx
        TermsZ <- object$Terms$termsYz
        TermsX.extra <- object$Terms$termsYx.extra
        TermsZ.extra <- object$Terms$termsYz.extra
        formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
        formYz <- Forms$formYz
        times <- Data$data[[timeVar]]
        GQsurv <- if (object$control$GQsurv == "GaussKronrod") 
            gaussKronrod() else gaussLegendre(object$control$GQsurv.k)
        wk <- GQsurv$wk
        sk <- GQsurv$sk
        K <- length(sk)
        anyLeftTrunc <- object$y$anyLeftTrunc
        if (anyLeftTrunc) {
            ni <- tapply(object$y$id, object$y$id, length)
            TimeL <- rep(object$y$TimeL, ni)
            P <- (times - TimeL) / 2
            st <- outer(P, sk) + c(times + TimeL) / 2
        } else {
            P <- times / 2
            st <- outer(P, sk + 1)
        }
        id.GK <- rep(seq_along(times), each = K)
        indBetas <- object$y$indBetas
        data.id2 <- Data$data.id[rep(object$y$id, each = K), ]
        data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
        if (param %in% c("td-value", "td-both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            Xs <- model.matrix(formYx, mfX)
            Zs <- model.matrix(formYz, mfZ)
        }
        if (param %in% c("td-extra", "td-both")) {
            mfX.extra <- model.frame(TermsX.extra, data = data.id2)
            mfZ.extra<- model.frame(TermsZ.extra, data = data.id2)
            Xs.extra <- model.matrix(Forms$extraForm$fixed, mfX.extra)
            Zs.extra <- model.matrix(Forms$extraForm$random, mfZ.extra)
        }
        betas <- object$postMeans$betas
        sigma <- object$postMeans$sigma
        D <- object$postMeans$D
        gammas <- object$postMeans$gammas
        alphas <- object$postMeans$alphas
        Dalphas <- object$postMeans$Dalphas
        if (nullY) {
            alphas <- rep(0, length.out = length(alphas))
            Dalphas <- rep(0, length.out = length(Dalphas))
        }
        Bs.gammas <- object$postMeans$Bs.gammas
        b <- ranef(object)
        idK <- rep(object$y$id, each = K)
        b <- b[idK, ]
        if (param %in% c("td-value", "td-both")) {
            Ys <- Funs$transFun.value(as.vector(Xs %*% betas + rowSums(Zs * b)), data.id2)
        }
        if (param %in% c("td-extra", "td-both")) {
            Ys.extra <- Funs$transFun.extra(as.vector(Xs.extra %*% betas[indFixed]) + 
                rowSums(Zs.extra * b[, indRandom, drop = FALSE]), data.id2)
        }
        tt <- c(switch(param,
                     "td-value" = as.matrix(Ys) %*% alphas, 
                     "td-extra" = as.matrix(Ys.extra) %*% Dalphas,
                     "td-both" = as.matrix(Ys) %*% alphas +  as.matrix(Ys.extra) %*% Dalphas,
                     "shared-betasRE" = (rep(betas[indBetas], each = nrow(b)) + b) %*% alphas,
                     "shared-RE" = b %*% alphas))
        W <- object$x$W[object$y$id, seq_along(gammas), drop = FALSE]
        eta.tw <- if (!is.null(W)) c(W %*% gammas) else rep(0, length(object$y$id))
        kn <- object$control$knots
        W2s <- splineDesign(unlist(kn, use.names = FALSE), c(t(st)), 
                            ord = object$control$ordSpline, outer.ok = TRUE)
        Vi <- exp(c(W2s %*% Bs.gammas) + tt)
        cumHaz <- exp(eta.tw) * P * tapply(rep(wk, length.out = length(Vi)) * Vi, id.GK, sum)
        names(cumHaz) <- row.names(Data$data)
        cumHaz
    }
}
