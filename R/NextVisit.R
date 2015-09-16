NextVisit <-
function (object, newdata, Tstart, Dt, tht, by = Dt/2, K = 10,
                       idVar = "id", onlyData = FALSE, betas = NULL, 
                       sigma = NULL, b = NULL, fun = NULL, simulate = FALSE, M = 300, 
                       seed = 123L) {
    if (!inherits(object, "JMbayes"))
        stop("'object' must inherit from class JMbayes.")
    timeVar <- object$timeVar
    times <- newdata[[timeVar]]
    TermsT <- object$Terms$termsT
    TimeVar <- all.vars(TermsT)[1L]
    Time <- newdata[[TimeVar]]
    timeSeq <- Tstart + 0.5 * Dt
    ntimeSeq <- length(timeSeq)
    newdata2 <- newdata[Time > Tstart & times <= Tstart, ]
    newdata2[[idVar]] <- newdata2[[idVar]][, drop = TRUE]
    set.seed(seed)
    sfit <- survfitJM(object, newdata2, survTimes = Tstart + Dt, 
                      simulate = simulate, M = M, seed = seed, idVar = idVar)
    pi.u.t <- sapply(sfit$summaries, "[", 2L)
    ind <- pi.u.t <= tht
    data.i <- split(newdata2, newdata2[[idVar]])
    data.i <- data.i2 <- data.i[sapply(data.i, nrow) > 0L]
    data.i <- data.i[ind]
    if (any(ind)) {
        newdata2Orig <- do.call(rbind, data.i)
        Time <- newdata2Orig[[TimeVar]]
        newdata2Orig <- newdata2Orig[Time > Tstart + 0.50001 * Dt, ]
        cvDCLOrig <- if (nrow(newdata2Orig)) {
            ss <- survfitJM(object, newdata = newdata2Orig, type = "Density",
                            idVar = idVar, simulate = simulate, M = M,
                            last.time = Tstart + 0.50001 * Dt, init.b = sfit$modes.b,
                            seed = seed)
            condDensOrig <- matrix(exp(-unlist(ss$full.results)), length(ss$summaries))
            condDensOrig[!is.finite(condDensOrig)] <- as.numeric(NA)
            condDensOrig <- apply(condDensOrig, 1, median, na.rm = TRUE)
            sum(-log(condDensOrig))
        } else as.numeric(NA)
        TermsX <- object$Terms$termsYx
        TermsZ <- object$Terms$termsYz
        formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
        formYz <- object$Forms$formYz
        if (is.null(betas))
            betas <- object$mcmc$betas
        if (is.null(sigma))
            sigma <- object$mcmc$sigma
        if (is.null(b))
            b <- sfit$modes.b[ind, , drop = FALSE]
        respVar <- all.vars(TermsX)[1L]
        if (is.null(fun) || !is.function(fun))
            fun <- function (x) x
        rand <- sample(nrow(betas), min(K, nrow(betas)))
        nn <- length(unique(newdata2Orig[[idVar]]))
        condDens <- matrix(0, nn, K)
        for (k in seq_len(K)) {
            betas.k <- betas[rand[k], ]
            sigma.k <- sigma[rand[k], ]
            data.i. <- data.i2[ind] <- mapply(function (d, bi) {
                last.row <- tail(d, n = 1L)[rep(1L, ntimeSeq), ]
                last.row[[timeVar]] <- timeSeq
                d <- rbind(d, last.row)
                mfX <- model.frame.default(TermsX, data = d)
                mfZ <- model.frame.default(TermsZ, data = d)
                X <- model.matrix.default(formYx, mfX)
                Z <- model.matrix.default(formYz, mfZ)
                mu <- c(X %*% betas.k + Z %*% bi)
                ii <- tail(1:nrow(d), ntimeSeq)
                d[[respVar]][ii] <- fun(rnorm(ntimeSeq, mu[ii], sigma.k))
                d
            }, d = data.i, bi = split(b, row(b)), SIMPLIFY = FALSE)        
            newdata2 <- do.call(rbind, data.i.)
            Time <- newdata2[[TimeVar]]
            newdata2 <- newdata2[Time > Tstart + 0.50001 * Dt, ]
            if (onlyData) {
                newdata2 <- do.call(rbind, data.i2)
                Time <- newdata2[[TimeVar]]
                newdata2 <- newdata2[Time > Tstart + 0.50001 * Dt, ]
                return(newdata2)
            }
            condDens[, k] <- if (nrow(newdata2)) {
                ss <- survfitJM(object, newdata = newdata2, type = "Density",
                                idVar = idVar, simulate = simulate, M = M,
                                last.time = Tstart + 0.50001 * Dt, init.b = sfit$modes.b,
                                seed = k)
                condDens.k <- matrix(exp(-unlist(ss$full.results)), length(ss$summaries))
                condDens.k[!is.finite(condDens.k)] <- as.numeric(NA)
                apply(condDens.k, 1, median, na.rm = TRUE)
            } else rep(as.numeric(NA), nrow(condDens))
        }
        #cvDCLExtra <- sum(-log(rowMeans(condDens)))
        cvDCLExtra <- sum(-log(apply(condDens, 1, median)))
        c(cvDCLOrig = cvDCLOrig, cvDCLExtra = cvDCLExtra, n = nn)
    } else {
        if (onlyData) newdata2 else c(NA, NA, NA)
    }
}
