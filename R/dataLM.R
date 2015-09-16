dataLM <-
function (data, Tstart, idVar = "id", respVar = "y", timeVar = "time", evTimeVar = "Time", 
                    summary = c("value", "slope", "area"), tranfFun = function (x) x) {
    if (!is.data.frame(data) || nrow(data) == 0L)
        stop("'data' must be a data.frame with more than one rows.\n")
    if (is.null(data[[idVar]]))
        stop("'idVar' not in 'data'.\n")
    if (is.null(data[[respVar]]))
        stop("'respVar' not in 'data'.\n")
    if (is.null(data[[timeVar]]))
        stop("'timeVar' not in 'data'.\n")
    if (is.null(data[[evTimeVar]]))
        stop("'evTimeVar' not in 'data'.\n")
    summary <- match.arg(summary)
    time <- data[[timeVar]]
    Time <- data[[evTimeVar]]
    ND <- data[Time > Tstart & time <= Tstart, ]
    f <- factor(ND[[idVar]], unique(ND[[idVar]]))
    if (summary == "value") {    
        ND[tapply(row.names(ND), f, tail, 1), ]
    } else if (summary == "slope") {
        do.call(rbind, lapply(split(ND, f), function (d) {
            d <- tail(d, 2)
            d$slope <- if (nrow(d) == 1) 0 else diff(tranfFun(d[[respVar]])) / diff(d[[timeVar]])
            tail(d, 1)
        }))
    } else {
        do.call(rbind, lapply(split(ND, f), function (d) {
            if (d[[timeVar]][1] != 0)
                d[[timeVar]][1] <- 0
            y <- tranfFun(d[[respVar]])
            t <- c(d[[timeVar]], Tstart)
            d$area <- sum(diff(t) * y)
            tail(d, 1)
        }))
        
    }    
}
