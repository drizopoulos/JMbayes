blockUpdate <-
function (old, logPost, logProp, meanProp, sigmaProp, logPost.old, logProp.old, df = 4) {
    isMat <- is.matrix(old)
    n <- if (isMat) nrow(old) else 1L
    u <- runif(n)
    proposed <- rmvt(n, meanProp, sigmaProp, df = df) # change RE list of Sigmas and mus
    logPost.prop <- logPost(proposed)
    logProp.prop <- logProp(proposed)
    logRatio <- logPost.prop + logProp.old - logPost.old - logProp.prop
    a <- pmin(exp(logRatio), 1)
    ind <- u <= a
    new <- if (isMat) {
        old[ind, ] <- proposed[ind, ]
        old
    } else {
        if (ind) proposed else old
    }
    list(new = new, logPost.old = ifelse(ind, logPost.prop, logPost.old), 
         logProp.old = ifelse(ind, logProp.prop, logProp.old))
}
