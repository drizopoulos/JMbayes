slice.tau <-
function (logPost, current.tau, step) {
    logy <- logPost(current.tau) - rexp(1L)
    tt <- runif(1L, 0, step)
    L <- current.tau - tt
    R <- current.tau + step - tt
    while (L > 0 && logPost(L) > logy) {
        L <- L - step
    }
    while (logPost(R) > logy) {
        R <- R + step
    }
    L <- max(0, L)
    repeat {
        new.tau <- runif(1L, L, R)
        new.logPost <- logPost(new.tau)
        if (new.logPost >= logy)
            break
        if (new.tau < current.tau) {
            L <- new.tau
        } else {
            R <- new.tau
        }
    }
    new.tau
}
