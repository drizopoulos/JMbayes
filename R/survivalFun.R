survivalFun <-
function (object, times, id = NULL, CI = FALSE, nullY = FALSE, nullCovs = FALSE) {
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
    n <- length(object$y$Time)
    if (is.null(id))
        id <- seq_len(n)
    if (length(times) == 1L)
        times <- rep(times, length.out = n)
    GQsurv <- if (object$control$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(object$control$GQsurv.k)
    wk <- GQsurv$wk
    sk <- GQsurv$sk
    K <- length(sk)
    P <- times/2
    st <- outer(P, sk + 1)
    id.GK <- rep(seq_along(times), each = K)
    indBetas <- object$y$indBetas
    data.id2 <- Data$data.id[rep(seq_len(n), each = K), ]
    data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame.default(TermsX, data = data.id2)
        mfZ <- model.frame.default(TermsZ, data = data.id2)
        Xs <- model.matrix.default(formYx, mfX)
        Zs <- model.matrix.default(formYz, mfZ)
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame.default(TermsX.extra, data = data.id2)
        mfZ.extra<- model.frame.default(TermsZ.extra, data = data.id2)
        Xs.extra <- model.matrix.default(Forms$extraForm$fixed, mfX.extra)
        Zs.extra <- model.matrix.default(Forms$extraForm$random, mfZ.extra)
    }
    betas <- object$postMeans$betas
    sigma <- object$postMeans$sigma
    D <- object$postMeans$D
    gammas <- object$postMeans$gammas
    if (nullCovs) 
        gammas <- rep(0, length.out = length(gammas))
    alphas <- object$postMeans$alphas
    Dalphas <- object$postMeans$Dalphas
    if (nullY) {
        alphas <- rep(0, length.out = length(alphas))
        Dalphas <- rep(0, length.out = length(Dalphas))
    }
    Bs.gammas <- object$postMeans$Bs.gammas
    b <- ranef(object)
    idK <- rep(1:n, each = K)
    ff <- function (betas, sigma, D, gammas, alphas, Dalphas, Bs.gammas, b) {
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
        W <- object$x$W[, seq_along(gammas), drop = FALSE]
        eta.tw <- if (!is.null(W)) c(W %*% gammas) else rep(0, length(object$y$id))
        kn <- object$control$knots
        W2s <- splineDesign(unlist(kn, use.names = FALSE), c(t(st)), 
                            ord = object$control$ordSpline, outer.ok = TRUE)
        Vi <- exp(c(W2s %*% Bs.gammas) + tt)
        cumHaz <- exp(eta.tw) * P * tapply(rep(wk, length.out = length(Vi)) * Vi, id.GK, sum)
        exp(-cumHaz)
    }
    est <- ff(betas, sigma, D, gammas, alphas, Dalphas, Bs.gammas, b)[id]
    if (CI) {
        M <- nrow(object$mcmc$betas)
        MCMCests <- matrix(0, n, M)
        for (m in seq_len(M)) {
            betas <- object$mcmc$betas[m, ]
            sigma <- object$mcmc$sigma[m, ]
            D <- matrix(object$mcmc$D[m, ], ncol(b), ncol(b))
            gammas <- object$mcmc$gammas[m, ]
            if (nullCovs) 
                gammas <- rep(0, length.out = length(gammas))
            alphas <- object$mcmc$alphas[m, ]
            Dalphas <- object$mcmc$Dalphas[m, ]
            if (nullY) {
                alphas <- rep(0, length.out = length(alphas))
                Dalphas <- rep(0, length.out = length(Dalphas))
            }
            Bs.gammas <- object$mcmc$Bs.gammas[m, ]
            b <- object$mcmc$b[, , m]
            MCMCests[, m] <- ff(betas, sigma, D, gammas, alphas, Dalphas, Bs.gammas, b)
        }
        attr(est, "low") <- apply(MCMCests, 1L, quantile, probs = 0.025)[id]
        attr(est, "upp") <- apply(MCMCests, 1L, quantile, probs = 0.975)[id]
    }
    est
}
