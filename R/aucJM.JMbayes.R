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
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
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
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + 1e-07
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
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & di == 0) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
        ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- paste(names(Ti), names(Tj), sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        
        
        ##### Generate future predictions for ind2, 3, 4 where some events occur after dThoriz
        # Create indexes to organise and unpack predictions
        # Note unlist prevents sapply retunrning a list when any index is zero
        nams_pi2i <- unlist( sapply( strsplit(names(ind2[ind2]), "_"), "[", 1) )
        nams_pi3j <- unlist( sapply( strsplit(names(ind3[ind3]), "_"), "[", 2) )
        nams_pi4i <- unlist( sapply( strsplit(names(ind4[ind4]), "_"), "[", 1) )
        nams_pi4j <- unlist( sapply( strsplit(names(ind4[ind4]), "_"), "[", 2) )
        nams_to_pred <- unique( c(nams_pi2i, nams_pi3j, nams_pi4i, nams_pi4j))

        # Predictions conditional observed events
        if(length(nams_to_pred) > 0){
            cond_preds <- if (is_counting) {
                survfitJM(object, newdata = newdata2[id %in% nams_to_pred, ], idVar = idVar,
                            last.time = Time[nams_to_pred], survTimes = Thoriz,
                            simulate = simulate, M = M, LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = newdata2[id %in% nams_to_pred, ], idVar = idVar,
                          last.time = Time[nams_to_pred], survTimes = Thoriz,
                          simulate = simulate, M = M)
            }
        }
        
        if (any(ind2)) {
            # Extract ind2 predictions
            pi2 <- cond_preds$summaries[ match(nams_pi2i, names(cond_preds$summaries)) ]
            pi2 <- 1 - sapply(pi2, "[", 1, 2)
            ind[ind2] <- ind[ind2] * pi2
        }
        if (any(ind3)) {
            # Extract ind3 predictions
            pi3 <- cond_preds$summaries[ match(nams_pi3j, names(cond_preds$summaries)) ]
            pi3 <- sapply(pi3, "[", 1, 2)
            ind[ind3] <- ind[ind3] * pi3
        }
        if (any(ind4)) {
            # Extract ind4 predictions
            pi4_i <- cond_preds$summaries[ match(nams_pi4i, names(cond_preds$summaries)) ]
            pi4_i <- 1 - sapply(pi4_i, "[", 1, 2)
            pi4_j <- cond_preds$summaries[ match(nams_pi4j, names(cond_preds$summaries)) ]
            pi4_j <- sapply(pi4_j, "[", 1, 2)
            ind[ind4] <- ind[ind4] * pi4_i * pi4_j
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
