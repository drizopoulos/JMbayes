log.posterior.b <-
function (b, y, Mats, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    ids.i <- ids %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new + Z.i %*% b)
    logY <- densLong(y[id.i], mu.y, sigma.new, log = TRUE, data = newdata[id.i, ])
    log.p.yb <- sum(logY)
    log.p.b <- densRE(b, mu = rep(0, ncol(Z.i)), D = D.new, log = TRUE, prop = FALSE)
    MM <- Mats[[ii]]
    st <- MM$st
    st2 <- MM$st2
    wk <- MM$wk
    wk2 <- MM$wk2
    P <- MM$P
    P2 <- MM$P2
    W2s <- MM$W2s
    Xs <- MM$Xs
    Zs <- MM$Zs
    Xs.extra <- MM$Xs.extra
    Zs.extra <- MM$Zs.extra
    Xu <- MM$Xu
    Zu <- MM$Zu
    ind <- MM$ind
    idT <- MM$idT
    id.GK2 <- MM$id.GK2
    if (param %in% c("td-value", "td-both"))
        Ys <- transFun.value(c(Xs %*% betas.new + Zs %*% b), data.s[ids.i, ])
    if (param %in% c("td-extra", "td-both"))
        Ys.extra <- transFun.extra(c(Xs.extra %*% betas.new[indFixed] + 
                                         Zs.extra %*% b[indRandom]), data.s[ids.i, ])
    if (estimateWeightFun) {
        wFun <- wk2 * weightFun(st2, shapes.new, max.time)
        Yu <- transFun.value(P2 * fastSumID(wFun * c(Xu %*% betas.new + Zu %*% b), id.GK2), 
                             data.s[ids.i, ])
    }
    tt <- if (!estimateWeightFun) {
        c(switch(param,
                 "td-value" = as.matrix(Ys) %*% alphas.new, 
                 "td-extra" =  as.matrix(Ys.extra) %*% Dalphas.new,
                 "td-both" = as.matrix(Ys) %*% alphas.new + 
                       as.matrix(Ys.extra) %*% Dalphas.new,
                 "shared-betasRE" = rep(sum((betas[indBetas] + b) * alphas.new), length(st)),
                 "shared-RE" = rep(sum(b * alphas.new), length(st))))
    } else {
        c(as.matrix(Yu) %*% alphas.new)
    }
    eta.tw <- if (!is.null(W)) {
            as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    } else 0
    Vi <- exp(c(W2s %*% Bs.gammas.new) + tt)
    log.survival <- - sum(exp(eta.tw) * P * fastSumID(wk * Vi, idT))
    if (all(st == 0))
        log.survival <- 1
    log.p.yb + log.survival + log.p.b
}
