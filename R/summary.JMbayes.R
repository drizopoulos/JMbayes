summary.JMbayes <- function (object, include.baseHazCoefs = FALSE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    coefs <- object$postMeans
    strs <- object$StErr
    sds <- object$StDev
    CIs <- object$CIs
    Ps <- object$Pvalues
    coefsY <- cbind("Value" = coefs$betas, "Std.Err" = strs$betas,
                    "Std.Dev" = sds$betas, "2.5%" = CIs$betas[1, ],
                    "97.5%" = CIs$betas[2, ],
                    "P" = Ps$betas)
    gammas <- c(coefs$gammas, coefs$alphas, coefs$Dalphas, coefs$shapes, coefs$Bs.gammas, coefs$tauBs)
    strs.gammas <- c(strs$gammas, strs$alphas, strs$Dalphas, strs$shapes, strs$Bs.gammas, strs$tauBs)
    sds.gammas <- c(sds$gammas, sds$alphas, sds$Dalphas, sds$shapes, sds$Bs.gammas, sds$tauBs)
    cis.gammas <- rbind(c(CIs$gammas[1, ], if (is.matrix(CIs$alphas)) CIs$alphas[1, ] else CIs$alphas[1],
                          if (is.matrix(CIs$Dalphas)) CIs$Dalphas[1, ] else CIs$Dalphas[1], CIs$shapes[1, ],
                          CIs$Bs.gammas[1, ],
                          if (object$baseHaz == "P-splines") CIs$tauBs[1]),
                        c(CIs$gammas[2, ], if (is.matrix(CIs$alphas)) CIs$alphas[2, ] else CIs$alphas[2],
                          if (is.matrix(CIs$Dalphas)) CIs$Dalphas[2, ] else CIs$Dalphas[2], CIs$shapes[2, ], CIs$Bs.gammas[2, ],
                          if (object$baseHaz == "P-splines") CIs$tauBs[2]))
    p.gammas <- c(Ps$gammas, Ps$alphas, Ps$Dalphas, Ps$shapes, Ps$Bs.gammas, if (!is.null(Ps$tauBs)) NA else NULL)
    coefsT <- cbind("Value" = gammas, "Std.Err" = strs.gammas, "Std.Dev" = sds.gammas,
                    "2.5%" = cis.gammas[1, ], "97.5%" = cis.gammas[2, ], "P" = p.gammas)
    if (!include.baseHazCoefs) {
        rnams <- rownames(coefsT)
        coefsT <- coefsT[-grep("Bs.gammas", rnams, fixed = TRUE), , drop = FALSE]
    } 
    out <- list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT,
                D = coefs$D, sigma = coefs$sigma)
    out$N <- nrow(object$x$X)
    out$n <- length(object$y$Time)
    out$d <- object$y$event
    out$id <- object$id
    out$control <- object$control
    out$time <- object$time
    out$knots <- unique(object$knots)
    out$baseHaz <- object$baseHaz
    out$conv <- object$conv
    out$param <- object$param
    out$estimateWeightFun <- object$estimateWeightFun
    out$df.RE <- object$df.RE
    out$densLongCheck <- object$densLongCheck
    out$DIC <- object$DIC
    out$pD <- object$pD
    out$LPML <- object$LPML
    out$call <- object$call
    class(out) <- "summary.JMbayes"
    out
}
