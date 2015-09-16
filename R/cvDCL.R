cvDCL <-
function (object, newdata, Tstart, idVar = "id", M = 300L, seed = 123L) {
    if (!inherits(object, "JMbayes"))
        stop("'object' must inherit from class JMbayes.")
    timeVar <- object$timeVar
    times <- newdata[[timeVar]]
    TermsT <- object$Terms$termsT
    TimeVar <- all.vars(TermsT)[1L]
    Time <- newdata[[TimeVar]]
    newdata2 <- newdata[Time > Tstart & times <= Tstart, ]
    newdata2[[idVar]] <- newdata2[[idVar]][, drop = TRUE]
    log_dens <- survfitJM(object, newdata2, seed = seed, idVar = idVar, type = "Density",
                          M = M, log = TRUE)
    mat <- matrix(unlist(log_dens$full.results, use.names = FALSE), ncol = M)
    sum(- log(rowMedians(exp(- mat))))
}
