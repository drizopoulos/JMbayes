marglogLik2 <- function (thetas, Data, priors, temp = 1.0, fixed_tau_Bs_gammas = FALSE) {
    # Data
    event <- Data[["event"]]
    idGK_fast <- Data[["idGK_fast"]]
    W1 <- Data[["W1"]]
    W1s <- Data[["W1s"]]
    event_colSumsW1 <- Data[["event_colSumsW1"]]
    W2 <- Data[["W2"]]
    W2s <- Data[["W2s"]]
    event_colSumsW2 <- Data[["event_colSumsW2"]]
    Wlong <- Data[["Wlong"]]
    Wlongs <- Data[["Wlongs"]]
    event_colSumsWlong <- Data[["event_colSumsWlong"]]
    Pw <- Data[["Pw"]]
    # priors
    mean_Bs_gammas <- priors[["mean_Bs_gammas"]]
    Tau_Bs_gammas <- priors[["Tau_Bs_gammas"]]
    mean_gammas <- priors[["mean_gammas"]]
    Tau_gammas <- priors[["Tau_gammas"]]
    mean_alphas <- priors[["mean_alphas"]]
    Tau_alphas <- priors[["Tau_alphas"]]
    A_tau_Bs_gammas <- priors[["A_tau_Bs_gammas"]]
    B_tau_Bs_gammas <- priors[["B_tau_Bs_gammas"]]
    # parameters
    list_thetas <- thetas
    if (fixed_tau_Bs_gammas) {
        list_thetas <- list_thetas[!names(list_thetas) %in% "tau_Bs_gammas"]
        tau_Bs_gammas <- thetas$tau_Bs_gammas
    }
    if (!ncol(W2)) {
        list_thetas <- list_thetas[!names(list_thetas) %in% "gammas"]
    }
    vec_thetas <- unlist(as.relistable(list_thetas))
    fn <- function (thetas) {
        thetas <- relist(thetas, skeleton = list_thetas)
        Bs_gammas <- thetas[["Bs_gammas"]]
        gammas <- thetas[["gammas"]]
        alphas <- thetas[["alphas"]]
        if (!fixed_tau_Bs_gammas)
            tau_Bs_gammas <- exp(thetas[["tau_Bs_gammas"]])
        if (is.null(gammas)) {
            - logPosterior_nogammas(temp, event, idGK_fast, W1, W1s, Bs_gammas,
                           Wlong, Wlongs, alphas, Pw, mean_Bs_gammas, Tau_Bs_gammas,
                           mean_alphas, Tau_alphas, tau_Bs_gammas, A_tau_Bs_gammas, B_tau_Bs_gammas)
        } else {
            - logPosterior(temp, event, idGK_fast, W1, W1s, Bs_gammas, W2, W2s,
                           gammas, Wlong, Wlongs, alphas, Pw, mean_Bs_gammas, Tau_Bs_gammas,
                           mean_gammas, Tau_gammas, mean_alphas, Tau_alphas,
                           tau_Bs_gammas, A_tau_Bs_gammas, B_tau_Bs_gammas)
        }
    }
    gr <- function (thetas) {
        thetas <- relist(thetas, skeleton = list_thetas)
        Bs_gammas <- thetas[["Bs_gammas"]]
        gammas <- thetas[["gammas"]]
        alphas <- thetas[["alphas"]]
        if (!fixed_tau_Bs_gammas)
            tau_Bs_gammas <- exp(thetas[["tau_Bs_gammas"]])
        out <- if (is.null(gammas)) {
            gradient_logPosterior_nogammas(temp, event, idGK_fast, event_colSumsW1, W1s,
                                  Bs_gammas, event_colSumsWlong, Wlongs, alphas, Pw,
                                  mean_Bs_gammas, Tau_Bs_gammas, mean_alphas, Tau_alphas,
                                  tau_Bs_gammas, A_tau_Bs_gammas, B_tau_Bs_gammas)
        } else {
            gradient_logPosterior(temp, event, idGK_fast, event_colSumsW1, W1s,
                                  Bs_gammas, event_colSumsW2, W2s, gammas,
                                  event_colSumsWlong, Wlongs, alphas, Pw,
                                  mean_Bs_gammas, Tau_Bs_gammas, mean_gammas,
                                  Tau_gammas, mean_alphas, Tau_alphas, tau_Bs_gammas,
                                  A_tau_Bs_gammas, B_tau_Bs_gammas)
        }
        if (fixed_tau_Bs_gammas)
            out <- head(out, -1)
        - out
    }
    hes <- function (thetas) {
        cd.vec(thetas, gr)
    }
    cd.vec <- function (x, f, ..., eps = 0.001) {
        n <- length(x)
        res <- matrix(0, n, n)
        ex <- pmax(abs(x), 1)
        for (i in 1:n) {
            x1 <- x2 <- x
            x1[i] <- x[i] + eps * ex[i]
            x2[i] <- x[i] - eps * ex[i]
            diff.f <- c(f(x1, ...) - f(x2, ...))
            diff.x <- x1[i] - x2[i]
            res[, i] <- diff.f/diff.x
        }
        0.5 * (res + t(res))
    }
    d <- length(unlist(list_thetas))
    pscale <- if (fixed_tau_Bs_gammas) rep(1, d) else c(rep(1, d - 1), 0.1)
    opt <- optim(vec_thetas, fn, gr, method = "BFGS", hessian = TRUE,
                 control = list(parscale = pscale))
    log_det_hessian <- determinant(opt$hessian)$modulus
    out <- as.vector(0.5 * (d * log(2 * pi) - log_det_hessian) - opt$value)
    #iLap <- iLap(opt, ff = fn, ff.gr = gr, ff.hess = hes, control = list(n.cores = 6))
    if (fixed_tau_Bs_gammas) {
        invH <- solve(opt$hessian)
        ind <- rep(seq_along(list_thetas), sapply(list_thetas, length))
        if (ncol(W2)) {
            ind_Bs_gammas <- which(ind == 1)
            ind_gammas <- which(ind == 2)
            ind_alphas <- which(ind == 3)
            Covs <- list(Bs_gammas = invH[ind_Bs_gammas, ind_Bs_gammas, drop = FALSE],
                         gammas = invH[ind_gammas, ind_gammas, drop = FALSE],
                         alphas = invH[ind_alphas, ind_alphas, drop = FALSE])
        } else {
            ind_Bs_gammas <- which(ind == 1)
            ind_alphas <- which(ind == 2)
            Covs <- list(Bs_gammas = invH[ind_Bs_gammas, ind_Bs_gammas, drop = FALSE],
                         alphas = invH[ind_alphas, ind_alphas, drop = FALSE])
        }
        attr(out, "Covs") <- Covs
        attr(out, "inits") <- c(relist(opt$par, skeleton = list_thetas),
                                list("tau_Bs_gammas"= tau_Bs_gammas))
    }
    out
}
