survPriors <-
function (Time, event, W, long, b, param) {
    if (param == "time-dependent") {
        fit <- survreg(Surv(Time, event) ~ W[, -1] + long)
        coefs <- - fit$coefficients / fit$scale
        V <- vcov(fit)
        nam <- "long"
        ind <- !names(coefs) %in% nam
        list(gammas = coefs[ind], var.gammas = V[which(ind), which(ind)], 
             alphas = coefs[!ind], var.alphas = V[which(!ind), which(!ind)])
    } else if (param == "shared-RE") {
        fit <- survreg(Surv(Time, event) ~ W[, -1] + b)
        coefs <- - fit$coefficients / fit$scale
        V <- vcov(fit)
        nam <- paste("b", seq_len(ncol(b)), sep = "")
        ind <- !names(coefs) %in% nam
        list(gammas = coefs[ind], var.gammas = V[which(ind), which(ind)], 
             alphas = coefs[!ind], var.alphas = V[which(!ind), which(!ind)])
    } else {
        fit <- survreg(Surv(Time, event) ~ W[, -1])
        coefs <- - fit$coefficients / fit$scale
        V <- vcov(fit)
        ind <- seq_along(coefs)
        list(gammas = coefs, var.gammas = V[ind, ind])
    }
}
