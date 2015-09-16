adjustKRW <-
function (K, acceptRate, q) {
    low <- 0.4
    upp <- 0.7
    if (acceptRate >= low && acceptRate <= upp)
        return(K)
    if (acceptRate > upp)
        K <- max(q, K * 0.7)
    if (acceptRate < low)
        K <- max(q, K * 1.2)
    K
}
