adjustScaleRW <-
function (scale, acceptRate, d, batch = NULL, startScale = NULL) {
    if (d == 1L) {
        low <- 0.4
        upp <- 0.5        
    } else if (d > 1L && d < 5L) {
        low <- 0.2
        upp <- 0.3
    } else {
        low <- 0.15
        upp <- 0.3
    }
    #if (acceptRate >= low && acceptRate <= upp)
    #    return(scale)
    qq <- if (is.null(batch)) 0.3 else min(0.1, 1/sqrt(batch))
    if (acceptRate > upp)
        scale <- scale * (1 + qq)
    if (acceptRate < low)
        scale <- scale * (1 - qq)
    if (!is.null(startScale))
        scale <- pmax(pmin(scale, startScale*1.25), startScale*0.75)
    min(max(scale, 1e-04), 5.0)
}
