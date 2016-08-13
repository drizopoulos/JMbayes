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
    auc <- if (length(Time) > 1) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        pi.u.t.i <- pi.u.t[pairs[1, ]]
        pi.u.t.j <- pi.u.t[pairs[2, ]]
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti < Thoriz & di == 1) & (Tj == Thoriz & dj == 0)
        ind3 <- (Ti < Thoriz & di == 0) & Tj > Thoriz
        ind <- ind1 | ind2 | ind3
        if (any(ind3)) {
            nams <- unique(names(ind3[ind3]))
            newdata3 <- newdata2[id %in% nams, ]
            tt <- model.response(model.frame(TermsT, newdata3))[, 1]
            pi2 <- numeric(nrow(newdata3))
            for (l in seq_along(pi2)) {
                pi2[l] <- c(summary(survfit(object, newdata = newdata3[l, ]), times = Thoriz)$surv) /
                    c(summary(survfit(object, newdata = newdata3[l, ]), times = tt[l])$surv)
            }
            pi2 <- 1 - pi2
            ind[ind3] <- ind[ind3] * pi2
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
