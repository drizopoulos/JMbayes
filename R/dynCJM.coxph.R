dynCJM.coxph <-
function (object, newdata, Dt, idVar = "id", t.max = NULL, timeVar = "time", 
                          weightFun = NULL, respVar = "y", evTimeVar = "Time",
                          summary = c("value", "slope", "area"), tranfFun = function (x) x, ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (!is.numeric(Dt) && length(Dt) > 1)
        stop("'Dt' must be a numeric scalar.\n")
    if (!is.null(weightFun) && !is.function(weightFun))
        stop("'weightFun' must be a function.\n")
    newdata$area <- newdata$slope <- 0
    TermsT <- object$terms
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    event <- SurvT[, 2]
    if (is.null(t.max) || !is.numeric(t.max) || length(t.max) > 1)
        t.max <- max(Time) + 1e-05
    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    P <- t.max / 2
    st <- P * (sk + 1)
    k <- length(st)
    auc.st <- numeric(k)
    form <- as.formula(paste(as.character(formula(object))[c(2,1,3)], collapse = " "))
    old <- options(warn = (2))
    on.exit(options(old))
    for (i in 1:k) {
        tt <- try({
            #data.i <- newdata[Time > st[i] & newdata[[timeVar]] <= st[i], ]
            #f <- factor(data.i[[idVar]], unique(data.i[[idVar]]))
            #data.i <- data.i[tapply(row.names(data.i), f, tail, 1), ]
            data.i <- dataLM(newdata, Tstart = st[i], idVar, respVar, timeVar, evTimeVar,
                             summary, tranfFun)
            object.i <- coxph(form, data = data.i)
            aucJM(object.i, newdata = newdata, Tstart = st[i], 
                  Dt = Dt, timeVar = timeVar, idVar = idVar,
                  respVar = respVar, evTimeVar = evTimeVar,
                  summary = summary, tranfFun = tranfFun)$auc
        }, TRUE)
        auc.st[i] <- if (!inherits(tt, "try-error")) tt else NA
    }
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
