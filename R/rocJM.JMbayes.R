rocJM.JMbayes <-
function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, idVar = "id", 
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
    Time <- SurvT[, 1]
    timeVar <- object$timeVar
    newdata2 <- newdata[Time > Tstart, ]
    newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]    
    pi.u.t <- survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                        simulate = simulate, M = M)
    pi.u.t <- sapply(pi.u.t$summaries, "[", 1, 2)
    # extract event process information
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[!duplicated(id), 1]
    event <- SurvT[!duplicated(id), 2]
    names(Time) <- names(event) <- as.character(unique(id))
    # subjects who died before Thoriz
    ind1 <- Time <= Thoriz & event == 1
    # subjects who were censored in the interval (Tstart, Thoriz)
    ind2 <- Time <= Thoriz & event == 0
    ind <- ind1 | ind2
    if (any(ind2)) {
        nams <- unique(names(ind2[ind2]))
        pi2 <- survfitJM(object, newdata = newdata2[id %in% nams, ], idVar = idVar, 
                         last.time = Time[nams], survTimes = Thoriz, 
                         simulate = simulate, M = M)
        pi2 <- 1 - sapply(pi2$summaries, "[", 1, 2)
        nams2 <- names(ind2[ind2])
        ind[ind2] <- ind[ind2] * pi2[nams2]
    }
    # calculate sensitivity and specificity
    thrs <- seq(0, 1, length = 101)
    TP <- colSums(outer(pi.u.t, thrs, "<") * ind) / sum(ind)
    FP <- colSums(outer(pi.u.t, thrs, "<") * (1 - ind)) / sum(1 - ind)
    Q <- colMeans(outer(pi.u.t, thrs, "<"))
    Q. <- 1 - Q
    k.1.0 <- (TP - Q) / Q.
    k.0.0 <- (1 - FP - Q.) / Q
    P <- mean(ind)
    P. <- 1 - P
    k.05.0 <- (P * Q. * k.1.0 + P. * Q * k.0.0) / (P * Q. + P. * Q)
    out <- list(TP = TP, FP = FP, qSN = k.1.0, qSP = k.0.0, qOverall = k.05.0, 
                thrs = thrs, 
                Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)), 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "rocJM"
    out
}
