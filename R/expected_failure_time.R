expectedCondFailureTime <- function (object, newdata, idVar = "id", last.time = NULL, 
                                     maxPossibleFailureTime = NULL) {
    if (!inherits(object, "JMbayes") && !inherits(object, "mvJMbayes"))
        stop("Use only with 'JMbayes' or 'mvJMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0L)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(maxPossibleFailureTime)) {
        Time <- if (inherits(object, "JMbayes")) object$y$Time 
        else object$model_info$coxph_components$Time
        maxPossibleFailureTime <- max(Time) * 1.5
    }
    timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
    cond_survival <- function (futureTimes, object, newdata, last.time, idVar) {
        sfit <- survfitJM(object, newdata, last.time = last.time, idVar = idVar, 
                          survTimes = futureTimes)
        sfit$summaries[[1]][, "Mean"]
    }
    newdata_split <- split(newdata, factor(newdata[[idVar]]))
    n <- length(newdata_split)
    means <- numeric(n)
    for (i in seq_len(n)) {
        ND <- newdata_split[[i]]
        if (is.null(last.time)) {
            last.time <- max(ND[[timeVar]])
        }
        Int <- integrate(cond_survival, lower = last.time, object = object, 
                         idVar = idVar, upper = maxPossibleFailureTime, newdata = ND, 
                         last.time = last.time, rel.tol = 0.05)$value
        means[i] <- last.time + Int
    }
    means
}

