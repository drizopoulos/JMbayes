myt <- function(times = 1) {
    tb <- "    "
    paste(rep(tb, times), collapse = "")
}

outcome_define <- function (family, outcome) {
    dens <- switch(family$family,
                   "gaussian" = "dnorm(mu",
                   "binomial" = "dbin(mu",
                   "poisson" = "dpois(mu")
    i <- outcome
    ind <- function (i) paste0("[j", i, "]")
    out <- paste0("y", i, ind(i), " ~ ", dens, i, ind(i))
    if (family$family == "gaussian") {
       out <- paste0(out, ", tau", i, ")")
    } else if (family$family == "binomial") {
        out <- paste0(out, ", 1", ")")
    } else {
        out <- paste0(out, ")")
    }
    out
}

linear_predictor <- function (family, outcome, colmns_HC, colmns_nHC,
                              overdispersion = FALSE, n_RE) {
    if (!family$link %in% c("identity", "log", "logit", "probit", "cloglog")) {
        stop("not permissible link function.")
    }
    i <- outcome
    ind <- function (i) paste0("[j", i, "]")
    ind. <- function (i) paste0("[j", i, ", ]")
    if (length(colmns_HC) > 1) {
        right <- if (n_RE == 1) {
            paste0(" <- u[id", i, ind(i), "] *", "Z", i, ind.(i))
        } else {
            paste0(" <- inprod(u[id", i, ind(i), ", RE_ind", i, "], ",
                   "Z", i, ind.(i), ")")
        }
    } else {
        right <- if (n_RE == 1) {
            paste0(" <- u[id", i, ind(i), "] * ", "Z", i, ind(i))
        } else {
            paste0(" <- u[id", i, ind(i), ", RE_ind", i, "] * ", "Z", i, ind(i))
        }
    }
    if (length(colmns_nHC)) {
        jj <- colmns_nHC
        if (length(jj) > 1) {
            kk <- paste0("c(", paste(jj, collapse = ", "), ")")
            right <- paste(right, paste0("+ inprod(betas", i, "[", kk, "], ",
                                         "X", i, "[j", i, ", ", kk, "])"))

        } else {
            right <- paste(right, paste0("+ betas", i, "[", jj, "] * ",
                                         "X", i, "[j", i, ", ", jj, "]"))
        }
    }
    if (family$family == "poisson" && overdispersion) {
        right <- paste0(right, " + uu", i, ind(i))
        right <- paste0(right, "\n", myt(2), "uu", i, ind(i), " ~ dnorm(0.0, tau_uu", i, ")")
    }
    out <- if (family$link == "identity") {
        paste0("mu", i, ind(i), right)
    } else {
        paste0(paste0("eta", i, ind(i), right, "\n"),
               myt(2), family$link, "(mu", i, ind(i), ") <- max(-15, min(15, eta", i, ind(i), "))")
    }
    out
}

build_outcome <- function (family, outcome, colmns_HC, colmns_nHC, overdispersion, n_RE) {
    i <- outcome
    paste0(myt(), "for (j", i, " in 1:N", i, ") {\n",
           myt(2), linear_predictor(family, outcome, colmns_HC, colmns_nHC,
                                    overdispersion, n_RE), "\n",
           myt(2), outcome_define(family, outcome), "\n",
           myt(), "}\n")
}

hc_linpred <- function (outcome, colmns_HC, incr, n_RE) {
    out <- ""
    i <- outcome
    if (n_RE == 1) {
        jj <- colmns_HC[[1]]
        if (length(jj) > 1) {
            ind.cl <- paste0("c(", paste(jj, collapse = ", "), ")")
            out <- paste0(out, paste0(myt(2), "mu.u[i] <- inprod(betas",
                                      i, "[", ind.cl, "], Xhc", i, "[i, ", ind.cl, "])\n"))
        } else {
            out <- paste0(out, paste0(myt(2), "mu.u[i] <- betas",
                                      i, "[", jj, "]\n"))
        }
        out <- paste0(out, paste0(myt(2), "b[i] <- u[i] - mu.u[i]\n"))
    } else {
        for (j in seq_along(colmns_HC)) {
            jj <- colmns_HC[[j]]
            j.incr <- j + incr
            if (length(jj) > 1) {
                ind.cl <- paste0("c(", paste(jj, collapse = ", "), ")")
                out <- paste0(out, paste0(myt(2), "mu.u[i, ", j.incr, "] <- inprod(betas",
                                          i, "[", ind.cl, "], Xhc", i, "[i, ", ind.cl, "])\n"))
            } else {
                out <- paste0(out, paste0(myt(2), "mu.u[i, ", j.incr, "] <- betas",
                                          i, "[", jj, "]\n"))
            }
            out <- paste0(out, paste0(myt(2), "b[i, ", j.incr, "] <- u[i, ",
                                      j.incr, "] - mu.u[i, ", j.incr, "]\n"))
        }
    }
    out
}

build_random_effects <- function (outcomes, colmns_HC, n_RE) {
    ns <- cumsum(c(0, head(sapply(colmns_HC, length), -1)))
    re_spec <- if (n_RE == 1) {
        "u[i] ~ dnorm(mu.u[i], inv_D)\n"
    } else {
        "u[i, 1:n_RE] ~ dmnorm(mu.u[i, ], inv_D[, ])\n"
    }
    paste0(myt(), "for (i in 1:n) {\n",
           paste(mapply(hc_linpred, outcomes, colmns_HC, ns, MoreArgs = list(n_RE = n_RE)),
                 collapse = ""),
           myt(2), re_spec, myt(), "}\n")
}

priors_betas <- function (outcome) {
    i <- outcome
    paste0(myt(), "for (k", i, " in 1:ncx", i, ") {\n",
           myt(2), "betas", i, "[k", i, "] ~ dnorm(0.0, tau_betas", i, ")\n",
           myt(), "}\n")
}

priors_taus <- function (family, outcome, overdispersion) {
    i <- outcome
    if (family$family == "gaussian") {
        paste0(myt(), "tau", i, " ~ dgamma(A_tau, B_tau)\n",
               myt(), "sigma", i, " <- 1 / sqrt(tau", i, ")\n")
    } else if (family$family == "poisson" && overdispersion) {
        paste0(myt(), "tau_uu", i, " ~ dgamma(A_tau, B_tau)\n")
    } else ""
}

priors_VC <- function (colmns_HC) {
    if (length(unlist(colmns_HC)) == 1) {
        "   inv_D ~ dt(0.0, tau_half_cauchy, 1)T(0, )\n"
    } else {
        paste("   inv_D ~ dwish(4 * priorR_D[, ], priorK_D)\n",
              "    for (l in 1:n_RE) {\n",
              "        priorR_D[l, l] ~ dgamma(A_RD, B_RD)\n",
              "    }\n")
    }
}

build_model <- function (families, outcomes, colmns_HC, colmns_nHC, overdispersion, n_RE) {
    models <- paste(mapply(build_outcome, families, outcomes, colmns_HC, colmns_nHC,
                           MoreArgs = list(overdispersion = overdispersion, n_RE = n_RE)),
                    collapse = "")
    REs <- build_random_effects(outcomes, colmns_HC, n_RE)
    priors <- paste(mapply(priors_betas, outcomes), collapse = "")
    priors <- paste0(priors, paste(mapply(priors_taus, families, outcomes,
                                          MoreArgs = list(overdispersion = overdispersion)), collapse = ""))
    #priors <- paste(priors, "   inv.D ~ dwish(4 * priorR.D[, ], priorK.D)\n")
    #priors <- paste0(priors, "    for (l in 1:n_RE) {\n",
    #                 "        priorR.D[l, l] ~ dgamma(A_R.D, B_R.D)\n",
    #                 "    }\n")
    priors <- paste(priors, priors_VC(colmns_HC))
    paste0("model {\n", models, REs, priors, "}\n", collapse = "")
}






