find_thresholds <- function (object, newdata, Dt, ...) {
    UseMethod("find_thresholds")
}

find_thresholds.mvJMbayes <- function (object, newdata, Dt, idVar = "id", M = 200L, 
                             n_cores =  max(1, parallel::detectCores() - 2), ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0L)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    Time <- object$model_info$coxph_components$Time
    event <- object$model_info$coxph_components$event
    nevents <- sum(event)
    ss <- seq(0, 1, length.out = floor(nevents / 20) + 2)
    times <- quantile(Time[event == 1], probs = tail(head(ss, -1), -1))
    do_roc <- function (i, object, newdata, times, Dt, idVar, M) {
        roc <- rocJM(object, newdata = newdata, Tstart = times[i], Dt = Dt, 
                     idVar = idVar, M = M)
        c("F1score" = roc$F1score, "Youden" = roc$Youden)
    }
    registerDoParallel(n_cores)
    out <- foreach(i = seq_along(times), .packages = "JMbayes", .combine = rbind) %dopar% {
        do_roc(i, object, newdata, times, Dt, idVar, M)
    }
    stopImplicitCluster()
    out <- cbind(times, out)
    colnames(out) <- c("times", "F1score", "Youden")
    rownames(out) <- NULL
    class(out) <- "ROC_cutoff"
    out
}

predict_eventTime <- function (object, newdata, cutpoints, ...) {
    UseMethod("predict_eventTime")
}

predict_eventTime.mvJMbayes <- function (object, newdata, cut_points, idVar = "id", 
                                         M = 500L, low_percentile = 0.025, ...) {
    if (!inherits(object, "mvJMbayes"))
        stop("Use only with 'mvJMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0L)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    if (length(unique(newdata[[idVar]])) > 1)
        stop("'predict_eventTime()' currently works for a single subject in 'newdata'.\n")
    Time <- object$model_info$coxph_components$Time
    max_Time <- max(Time)
    last_time <- max(newdata[[object$model_info$timeVar]])
    sfit <- survfitJM(object, newdata = ND, M = M, idVar = idVar,
                      survTimes = seq(last_time, max_Time, length.out = 225))
    sfit <- sfit$summaries[[1]][, c('times', 'Mean')]
    extract_time <- function (sfit, percentile) {
        sfit[which.min(abs(sfit[, 2] - percentile)), 1]
    }
    low_time <- extract_time(sfit, 1 - low_percentile)
    median_time <- extract_time(sfit, 0.5)
    if (median_time - low_time < quantile(Time, 0.25)) {
        median_time
    } else {
        cut_last_time <- cut_points[which.min(abs(cut_points[, 1] - last_time)), 2]
        extract_time(sfit, cut_last_time)
    }
}