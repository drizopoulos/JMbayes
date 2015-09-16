slice.shape <-
function (logPost, current.shapes, step, which, iters = 6L) {
    current.shape <- current.shapes[which]
    logy <- logPost(current.shape, which)[[1L]] - rexp(1L)
    tt <- runif(1L, 0, step)
    L <- current.shape - tt
    R <- current.shape + step - tt
    count <- 0L
    lP.val <- logPost(L, which)[[1L]]
    while (L > 0 && !is.na(lP.val) && lP.val > logy) {
        L <- L - step
        count <- count + 1L
        if (count > iters)
            break
    }
    count <- 0L
    lP.val <- logPost(R, which)[[1L]]
    while (!is.na(lP.val) && lP.val > logy) {
        R <- R + step
        count <- count + 1L
        if (count > iters)
            break
    }
    L <- max(0, L, na.rm = TRUE)
    count <- 0L
    repeat {
        count <- count + 1L
        new.shape <- runif(1L, L, R)
        new.logPost <- logPost(new.shape, which)
        new.logPost.val <- new.logPost[[1L]]
        if (is.na(new.logPost.val) || new.logPost.val >= logy || count > iters)
            break
        if (new.shape < current.shape) {
            L <- new.shape
        } else {
            R <- new.shape
        }
    }
    c(list(new.shape = new.shape, fail = count > iters || is.na(new.logPost.val)), 
      new.logPost)
}
