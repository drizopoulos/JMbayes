aucJM.JMbayes <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, idVar = "id", 
        simulate = FALSE, M = 100, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$Terms$termsT
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    is_counting <- attr(SurvT, "type") == "counting"
    Time <- if (is_counting) {
        ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
    } else {
        SurvT[, 1]
    }
    timeVar <- object$timeVar
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- newdata2[Time[ordTime] > Tstart, ]
    newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]
    pi.u.t <- if (is_counting) {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                  simulate = simulate, M = M, LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                  simulate = simulate, M = M)
    }
    pi.u.t <- sapply(pi.u.t$summaries, "[", 1, 2)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    if (is_counting) {
        f <- factor(id, levels = unique(id))
        Time <- tapply(SurvT[, 2], f, tail, 1)
        event <- tapply(SurvT[, 3], f, tail, 1)
    } else{
        Time <- SurvT[!duplicated(id), 1]
        event <- SurvT[!duplicated(id), 2]
    }
    names(Time) <- names(event) <- as.character(unique(id))
    if (!all(names(pi.u.t) == names(Time)))
        stop("mismatch between 'Time' variable names and survival probabilities names.")
    auc <- if (length(Time) > 1) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        pi.u.t.i <- pi.u.t[pairs[1, ]]
        pi.u.t.j <- pi.u.t[pairs[2, ]]
        ind1 <- (Ti < Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <  Thoriz & di == 1) & (Tj == Thoriz & dj == 0)
        ind3 <- (Ti < Thoriz & di == 0) & Tj > Thoriz
        ind <- ind1 | ind2 | ind3
        if (any(ind3)) {
            nams <- unique(names(ind3[ind3]))
            pi2 <- if (is_counting) {
                survfitJM(object, newdata = newdata2[id %in% nams, ], idVar = idVar, 
                             last.time = Time[nams], survTimes = Thoriz, 
                             simulate = simulate, M = M, LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = newdata2[id %in% nams, ], idVar = idVar, 
                          last.time = Time[nams], survTimes = Thoriz, 
                          simulate = simulate, M = M)
            }
            pi2 <- 1 - sapply(pi2$summaries, "[", 1, 2)
            nams2 <- names(ind3[ind3])
            ind[ind3] <- ind[ind3] * pi2[nams2]
        }
        sum((pi.u.t.i < pi.u.t.j) * c(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE)
    } else {
        NA
    }
    out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)), 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "aucJM"
    out
}
