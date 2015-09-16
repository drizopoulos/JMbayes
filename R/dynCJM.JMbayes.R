dynCJM.JMbayes <-
function (object, newdata, Dt, idVar = "id", t.max = NULL, simulate = FALSE, M = 100, 
                            weightFun = NULL, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (!is.numeric(Dt) && length(Dt) > 1)
        stop("'Dt' must be a numeric scalar.\n")
    if (!is.null(weightFun) && !is.function(weightFun))
        stop("'weightFun' must be a function.\n")
    TermsT <- object$Terms$termsT
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    event <- SurvT[, 2]
    if (is.null(t.max) || !is.numeric(t.max) || length(t.max) > 1)
        t.max <- max(Time) + 1e-05
    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    P <- t.max / 2
    st <- P * (sk + 1)
    auc.st <- sapply(st, function (t) {
        test <- try(aucJM(object, newdata = newdata, Tstart = t, Dt = Dt, idVar = idVar, simulate = simulate, M = M)$auc, TRUE)
        if (!inherits(test, "try-error")) test else NA
        
    })
    if (is.null(weightFun)) {
        weightFun <- function (t, Dt) {
            sfit <- survfit(Surv(Time, event) ~ 1)
            S.t <- summary(sfit, times = t)$surv
            S.tdt <- summary(sfit, times = t + Dt)$surv
            r <- (S.t - S.tdt) * S.tdt
            if (length(r)) r else NA
        }
    }
    w.st <- sapply(st, function (t) weightFun(t, Dt))
    dynC <- sum(wk * auc.st * w.st, na.rm = TRUE) / sum(wk * w.st, na.rm = TRUE)
    out <- list(dynC = dynC, times = st, AUCs = auc.st, weights = w.st, t.max = t.max, Dt = Dt, 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "dynCJM"
    out
}
