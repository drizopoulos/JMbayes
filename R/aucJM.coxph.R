aucJM.coxph <-
function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
          idVar = "id", respVar = "y", timeVar = "time", evTimeVar = "Time",
          summary = c("value", "slope", "area"), tranfFun = function (x) x, ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
    newdata$area <- newdata$slope <- 0
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$terms
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- dataLM(newdata2, Tstart, idVar, respVar, timeVar, evTimeVar, summary, 
                      tranfFun)
    pi.u.t <- c(summary(survfit(object, newdata = newdata2), times = Thoriz)$surv)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[!duplicated(id), 1]
    event <- SurvT[!duplicated(id), 2]
    names(pi.u.t) <- names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + 1e-07
    }
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
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            ND <- newdata2[id %in% unq_nams_i, ]
            tt <- model.response(model.frame(TermsT, ND))[, 1]
            pi2 <- numeric(nrow(ND))
            for (l in seq_along(pi2)) {
                obj <- survfit(object, newdata = ND[l, ])
                pi2[l] <- summary(obj, times = Thoriz)$surv / summary(obj, times = tt[l])$surv
            }
            pi2 <- 1 - pi2
            names(pi2) <- as.character(ND[[idVar]])
            ind[ind2] <- ind[ind2] * pi2[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            ND <- newdata2[id %in% unq_nams_j, ]
            tt <- model.response(model.frame(TermsT, ND))[, 1]
            pi3 <- numeric(nrow(ND))
            for (l in seq_along(pi3)) {
                obj <- survfit(object, newdata = ND[l, ])
                pi3[l] <- summary(obj, times = Thoriz)$surv / summary(obj, times = tt[l])$surv
            }
            names(pi3) <- as.character(ND[[idVar]])
            ind[ind3] <- ind[ind3] * pi3[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            ND_i <- newdata2[id %in% unq_nams_i, ]
            ND_j <- newdata2[id %in% unq_nams_j, ]
            tt_i <- model.response(model.frame(TermsT, ND_i))[, 1]
            tt_j <- model.response(model.frame(TermsT, ND_j))[, 1]
            pi4_i <- numeric(nrow(ND_i))
            for (l in seq_along(pi4_i)) {
                obj <- survfit(object, newdata = ND_i[l, ])
                pi4_i[l] <- summary(obj, times = Thoriz)$surv / summary(obj, times = tt_i[l])$surv
            }
            pi4_i <- 1 - pi4_i
            names(pi4_i) <- as.character(ND_i[[idVar]])
            pi4_j <- numeric(nrow(ND_j))
            for (l in seq_along(pi4_j)) {
                obj <- survfit(object, newdata = ND_j[l, ])
                pi4_j[l] <- summary(obj, times = Thoriz)$surv / summary(obj, times = tt_j[l])$surv
            }
            names(pi4_j) <- as.character(ND_j[[idVar]])
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
