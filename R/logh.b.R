logh.b <-
function (b, mats) {
    if (!is.matrix(b))
        b <- rbind(b)
    W2 <- mats$W2
    Xtime <- mats$Xtime
    Ztime <- mats$Ztime
    Xtime.extra <- mats$Xtime.extra
    Ztime.extra <- mats$Ztime.extra
    Xu <- mats$Xu
    Zu <- mats$Zu
    P <- mats$P
    st <- mats$st
    wk <- mats$wk
    id.GK <- mats$id.GK
    data <- mats$data
    logh0 <- as.vector(W2 %*% Bs.gammas.new)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas.new) else rep(0, length(logh0))
    if (param %in% c("td-value", "td-both") && !estimateWeightFun)
        Y <- transFun.value(c(Xtime %*% betas.new) + rowSums(Ztime * b), data)
    if (param %in% c("td-extra", "td-both"))
        Y.extra <- transFun.extra(c(Xtime.extra %*% betas.new[indFixed]) + 
                                    rowSums(Ztime.extra * b[, indRandom, drop = FALSE]), 
                                  data)
    if (estimateWeightFun) {
        wFun <- wk * weightFun(st, shapes.new, max.time)
        Yu <- transFun.value(P * fastSumID(wFun * (c(Xu %*% betas.new) + 
                            rowSums(Zu * b[id.GK, , drop = FALSE])), id.GK), data)
    }
    tt <- if (!estimateWeightFun) {
        c(switch(param, 
                 "td-value" = as.matrix(Y) %*% alphas.new, 
                 "td-extra" =  as.matrix(Y.extra) %*% Dalphas.new,
                 "td-both" = as.matrix(Y) %*% alphas.new + as.matrix(Y.extra) %*% Dalphas.new,
                 "shared-betasRE" = (b + rep(betas, each = nrow(b))) %*% alphas.new,
                 "shared-RE" = b %*% alphas.new))
    } else {
        c(as.matrix(Yu) %*% alphas.new)
    }
    logh0 + eta.tw + tt
}
