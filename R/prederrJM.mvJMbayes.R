prederrJM.mvJMbayes <- function (object, newdata, Tstart, Thoriz, lossFun = c("square", "absolute"), 
                               interval = FALSE, idVar = "id", M = 100, ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    lossFun <- if (is.function(lossFun)) {
        lf <- lossFun
        match.fun(lossFun)
    } else {
        lf <- match.arg(lossFun)
        if (lf == "absolute") function (x) abs(x) else function (x) x*x
    }
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$model_info$coxph_components$Terms
    environment(TermsT) <- .GlobalEnv
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    is_counting <- attr(SurvT, "type") == "counting"
    Time <- if (is_counting) {
        ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
    } else {
        SurvT[, 1]
    }
    timeVar <- object$model_info$timeVar
    newdata2 <- newdata[Time > Tstart, ]
    SurvT <- model.response(model.frame(TermsT, newdata2))
    if (is_counting) {
        id2 <- newdata2[[idVar]]
        f <- factor(id2, levels = unique(id2))
        Time <- ave(SurvT[, 2], f, FUN = function (x) tail(x, 1))
        delta <- ave(SurvT[, 3], f, FUN = function (x) tail(x, 1))
    } else {
        Time <- SurvT[, 1]
        delta <- SurvT[, 2]
    }
    timesInd <- newdata2[[timeVar]] <= Tstart
    aliveThoriz <- newdata2[Time > Thoriz & timesInd, ]
    deadThoriz <- newdata2[Time <= Thoriz & (delta == 1 | delta == 3) & timesInd, ]
    indCens <- Time < Thoriz & delta == 0 & timesInd
    censThoriz <- newdata2[indCens, ]
    nr <- length(unique(newdata2[[idVar]]))
    idalive <- unique(aliveThoriz[[idVar]])
    iddead <- unique(deadThoriz[[idVar]])
    idcens <- unique(censThoriz[[idVar]])
    prederr <- if (length(unique(Time)) > 1 && nrow(aliveThoriz) > 1 &&
                   nrow(deadThoriz) > 1) {
        Surv.aliveThoriz <- if (is_counting) {
            survfitJM(object, newdata = aliveThoriz, idVar = idVar, M = M,
                      survTimes = Thoriz, last.time = rep(Tstart, length(idalive)),
                      LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
            survfitJM(object, newdata = aliveThoriz, idVar = idVar, M = M,
                      survTimes = Thoriz, last.time = rep(Tstart, length(idalive)))
            
        }
        Surv.deadThoriz <- if (is_counting) {
            survfitJM(object, newdata = deadThoriz, idVar = idVar, 
                      survTimes = Thoriz, last.time = rep(Tstart, length(iddead)),
                      LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
            survfitJM(object, newdata = deadThoriz, idVar = idVar, 
                      survTimes = Thoriz, last.time = rep(Tstart, length(iddead)))
        }
        Surv.aliveThoriz <- sapply(Surv.aliveThoriz$summaries, "[", 2)
        Surv.deadThoriz <- sapply(Surv.deadThoriz$summaries, "[", 2)
        if (nrow(censThoriz)) {
            Surv.censThoriz <- if (is_counting) {
                survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                          survTimes = Thoriz, last.time = rep(Tstart, length(idcens)),
                          LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                          survTimes = Thoriz, last.time = rep(Tstart, length(idcens)))
            }
            tt <- Time[indCens]
            weights <- if (is_counting) {
                survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                          survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])],
                          LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                          survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])])
            }
            Surv.censThoriz <- sapply(Surv.censThoriz$summaries, "[", 2)
            weights <- sapply(weights$summaries, "[", 2)
        } else {
            Surv.censThoriz <- weights <- NA
        }
        if (!interval) {
            (1/nr) * sum(lossFun(1 - Surv.aliveThoriz), lossFun(0 - Surv.deadThoriz),
                         weights * lossFun(1 - Surv.censThoriz) + (1 - weights) * lossFun(0 - Surv.censThoriz), 
                         na.rm = TRUE)
        } else {
            TimeCens <- object$model_info$coxph_components$Time
            deltaCens <- 1 - object$model_info$coxph_components$event
            KMcens <- survfit(Surv(TimeCens, deltaCens) ~ 1)
            times <- TimeCens[TimeCens > Tstart & TimeCens < Thoriz & !deltaCens]
            times <- sort(unique(times))
            k <- as.numeric(table(times))
            w <- summary(KMcens, times = Tstart)$surv / summary(KMcens, times = times)$surv
            prederr.times <- sapply(times, 
                                    function (t) prederrJM(object, newdata, Tstart, t, M = M,
                                                           interval = FALSE, idVar = idVar)$prederr)
            num <- sum(prederr.times * w * k, na.rm = TRUE)
            den <- sum(w * k, na.rm = TRUE)
            num / den
        }
    } else {
        nr <- NA
        NA
    }
    out <- list(prederr = prederr, nr = nr, Tstart = Tstart, Thoriz = Thoriz, 
                interval = interval, classObject = class(object), 
                nameObject = deparse(substitute(object)), lossFun = lf)
    class(out) <- "prederrJM"
    out
}
