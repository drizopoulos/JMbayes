S.b <-
function (t, b, ii, Mats, log = FALSE) {
    if (t == 0) {
        if (log) return(0) else return(0)
    }
    idT.i <- idT %in% ii
    ids.i <- ids %in% ii
    st <- Mats$st
    st2 <- Mats$st2
    wk <- Mats$wk
    wk2 <- Mats$wk2
    P <- Mats$P
    P2 <- Mats$P2
    Xs <- Mats$Xs
    Zs <- Mats$Zs
    Xs.extra <- Mats$Xs.extra
    Zs.extra <- Mats$Zs.extra
    Xu <- Mats$Xu
    Zu <- Mats$Zu
    W2s <- Mats$W2s
    ind <- Mats$ind
    idT <- Mats$idT
    id.GK2 <- Mats$id.GK2
    if (param %in% c("td-value", "td-both") && !estimateWeightFun)
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
                 "td-both" = as.matrix(Ys) %*% alphas.new + as.matrix(Ys.extra) %*% Dalphas.new,
                 "shared-betasRE" = rep(sum((betas[indBetas] + b) * alphas.new), length(st)),
                 "shared-RE" = rep(sum(b * alphas.new), length(st))))
    } else {
        c(as.matrix(Yu) %*% alphas.new)
    }
    eta.tw <- if (!is.null(W)) {
        as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    } else 0
    Vi <- exp(c(W2s %*% Bs.gammas.new) + tt)
    ind <- Mats$st < min(Mats$kn)
    wk[ind] <- 0
    log.survival <- - sum(exp(eta.tw) * P * fastSumID(wk * Vi, idT))
    if (log) log.survival else exp(log.survival)
}
