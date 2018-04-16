aucJM.mvJMbayes <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, 
                             idVar = "id", M = 100, ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$model_info$coxph_components$Terms
    environment(TermsT) <- .GlobalEnv
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    is_counting <- attr(SurvT, "type") == "counting"
    is_interval <- attr(SurvT, "type") == "interval"
    Time <- if (is_counting) {
        ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
    } else if (is_interval) {
        Time1 <- SurvT[, "time1"]
        Time2 <- SurvT[, "time2"]
        Time <- Time1
        Time[Time2 != 1] <- Time2[Time2 != 1]
        Time
    } else {
        SurvT[, 1]
    }
    timeVar <- object$model_info$timeVar
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- newdata2[Time[ordTime] > Tstart, ]
    newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]
    pi.u.t <- if (is_counting) {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, M = M, 
                  LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, M = M)
    }
    pi.u.t <- sapply(pi.u.t$summaries, "[", 1, 2)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    if (is_counting) {
        f <- factor(id, levels = unique(id))
        Time <- tapply(SurvT[, 2], f, tail, 1)
        event <- tapply(SurvT[, 3], f, tail, 1)
    } else if (is_interval) {
        Time1 <- SurvT[, "time1"]
        Time2 <- SurvT[, "time2"]
        Time <- Time1
        Time[Time2 != 1] <- Time2[Time2 != 1]
        Time <- Time[!duplicated(id)]
        event <- SurvT[!duplicated(id), "status"]
    } else {
        Time <- SurvT[!duplicated(id), 1]
        event <- SurvT[!duplicated(id), 2]
    }
    names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + runif(length(Time[dupl]), 1e-07, 1e-06)
    }
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
        ind1 <- (Ti <= Thoriz & (di == 1 | di == 3)) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & (di == 0 | di == 2)) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & (di == 1 | di == 3)) & (Tj <= Thoriz & (dj == 0 | dj == 2))
        ind4 <- (Ti <= Thoriz & (di == 0 | di == 2)) & (Tj <= Thoriz & (dj == 0 | dj == 2))
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- paste(names(Ti), names(Tj), sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            pi2 <- if (is_counting) {
                survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], idVar = idVar, 
                             last.time = Time[unq_nams_i], survTimes = Thoriz, 
                             M = M, LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], idVar = idVar, 
                          last.time = Time[unq_nams_i], survTimes = Thoriz, M = M)
            }
            pi2 <- 1 - sapply(pi2$summaries, "[", 1, 2)
            ind[ind2] <- ind[ind2] * pi2[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            pi3 <- if (is_counting) {
                survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], idVar = idVar, 
                          last.time = Time[unq_nams_j], survTimes = Thoriz, 
                          M = M, LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], idVar = idVar, 
                          last.time = Time[unq_nams_j], survTimes = Thoriz, M = M)
            }
            pi3 <- sapply(pi3$summaries, "[", 1, 2)
            ind[ind3] <- ind[ind3] * pi3[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            if (is_counting) {
                pi4_i <- survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], 
                                   idVar = idVar, last.time = Time[unq_nams_i], survTimes = Thoriz, 
                                   M = M, 
                                   LeftTrunc_var = all.vars(TermsT)[1L])
                pi4_j <- survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], 
                                   idVar = idVar, last.time = Time[unq_nams_j], survTimes = Thoriz, 
                                   M = M, 
                                   LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                pi4_i <- survfitJM(object, newdata = newdata2[id %in% unq_nams_i, ], 
                                   idVar = idVar, last.time = Time[unq_nams_i], survTimes = Thoriz, 
                                   M = M)
                pi4_j <- survfitJM(object, newdata = newdata2[id %in% unq_nams_j, ], 
                                   idVar = idVar, last.time = Time[unq_nams_j], survTimes = Thoriz, 
                                   M = M)
            }
            pi4_i <- 1 - sapply(pi4_i$summaries, "[", 1, 2)
            pi4_j <- sapply(pi4_j$summaries, "[", 1, 2)
            ind[ind4] <- ind[ind4] * pi4_i[nams_i] * pi4_j[nams_j]
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
