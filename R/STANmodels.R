data_part <- function (families, ncolsZs, extraXs, n_RE, colmns_HC, colmns_nHC) {
    outcomes <- seq_along(families)
    data_outcome <- function (outcome, family, ncolsZ, extraX = FALSE, colmns_HC, colmns_nHC) {
        type_long_outcome <- switch(family$family,
                                    "gaussian" = paste0("vector[N", outcome, "] y", outcome),
                                    "binomial" = paste0("int<lower=0, upper=1> y", outcome, "[N", outcome, "]"),
                                    "poisson" = paste0("int<lower=0> y", outcome, "[N", outcome, "]")
        )
        f <- function (colmns_HC_i, outcome) {
            if (length(colmns_HC_i) > 1) {
                out <- vector("character", length(colmns_HC_i))
                for (k in seq_along(colmns_HC_i)) {
                    ll <- length(colmns_HC_i[[k]])
                    out[k] <- paste0(myt(), "int colmns_HC", outcome, k,
                                     if (ll > 1) paste0("[", ll, "]"), ";\n")
                }
                paste(out, collapse = "")
            } else {
                ll <- length(colmns_HC_i)
                paste0(myt(), "int colmns_HC", outcome,
                       if (ll > 1) paste0("[", ll, "]"), ";\n")
            }
        }
        paste0(myt(), "int N", outcome, ";\n",
               myt(), "int ncx", outcome, ";\n",
               myt(), "int id", outcome, "[N", outcome, "];\n",
               myt(), "int RE_ind", outcome, if (ncolsZ > 1) paste0("[", ncolsZ, "]"), ";\n",
               f(colmns_HC, outcome),
               if (ll <- length(colmns_nHC)) 
                   paste0(myt(), "int colmns_nHC", outcome, 
                          if (ll > 1) paste0("[", ll, "]"), ";\n"),
               myt(), type_long_outcome, ";\n", 
               myt(), "matrix[N", outcome, ", ", ncolsZ, "] ", "Z", outcome, ";\n",
               myt(), "matrix[n, ncx", outcome, "] ", "Xhc", outcome, ";\n",
               if (extraX) paste0(myt(), "matrix[N", outcome, ", ", "ncx", outcome, "] ",
                                  "X", outcome, ";\n"))
    }
    paste0("data {\n", myt(), "int n;\n", myt(), "int n_RE;\n", 
           paste0(mapply(data_outcome, outcomes, families, ncolsZs, extraXs, colmns_HC, colmns_nHC), 
                  collapse = ""),
           paste0(sapply(outcomes, function (outcome) 
               paste0(myt(), "real<lower=0> scale_betas", outcome, ";\n")), collapse = ""),
           if (any(sapply(families, `[[`, 'family') == "gaussian")) 
               paste0(myt(), "real<lower=0> scale_sigmas;\n"),
           myt(), "real<lower=0> scale_diag_D;\n",
           if (n_RE > 1)
               paste0(myt(), "real<lower=0> lkj_shape;\n"),
           "}\n")
}

##########################################################################################

parameters <- function (families, n_RE) {
    outcomes <- seq_along(families)
    set_parms <- function (outcome, family) {
        paste0(myt(), "vector[ncx", outcome, "] betas", outcome, ";\n",
               if (family$family == "gaussian") paste0(myt(), 
                                                       "real<lower = 0> sigma", outcome, ";\n"))
    }
    paste0("\nparameters {\n",
           paste0(mapply(set_parms, outcomes, families), collapse = ""),
           myt(), "matrix[n, n_RE] u;\n",
           if (n_RE > 1)
               paste0(myt(), "vector<lower = 0>[n_RE] L_var_D;\n",
                      myt(), "cholesky_factor_corr[n_RE] L_corr_D;\n")
           else 
               paste0(myt(), "real<lower = 0> D;\n"),
           "}\n")
}

##########################################################################################

transformed_parameters <- function (families, colmns_HC, colmns_nHC, RE_inds) {
    outcomes <- seq_along(families)
    def_etas <- function (outcome) {
        paste0(myt(), "vector[N", outcome, "] eta", outcome, ";\n")
    }
    colmns_HC2 <- unlist(colmns_HC, recursive = FALSE)
    nams <- names(colmns_HC2)
    ncols <- sapply(colmns_HC2, length)
    HC_part <- function (outcome, columns, nam, ncol, i) {
        #if (ncol > 1) {
        #    paste0(myt(2), "mu_u[i, ", i,"] = dot_product(Xhc", outcome, 
        #           "[i, ", nam, "], betas", outcome, "[", nam, "]);\n")
        #} else if (ncol == 1) {
        #    paste0(myt(2), "mu_u[i, ", i,"] = Xhc", outcome, 
        #           "[i, ", columns, "] * betas", outcome, "[", columns, "];\n")
        #} else {
        #    paste0(myt(2), "mu_u[i, ", i,"] = 0.0;\n")
        #}
        if (ncol) {
            paste0(myt(2), "mu_u[i, ", i,"] = ", 
                   paste0("Xhc", outcome, "[i, ", columns, "] * betas", 
                          outcome, "[", columns, "]", collapse = " + "), ";\n")
        } else {
            paste0(myt(2), "mu_u[i, ", i,"] = 0.0;\n")
        }
    }
    linpred_part <- function (outcome, colmns_HC, colmns_nHC, RE_ind) {
        j_i <- paste0("j", outcome)
        paste0(myt(), "for (", j_i, " in 1:N", outcome, ") {\n",
               myt(2), "eta", outcome, "[", j_i, "] = ", 
               #if (length(colmns_HC) > 1)
                #   paste0("dot_product(Z", outcome, "[", j_i,
                #          ", ], u[id", outcome, "[", j_i, "], RE_ind", outcome, "])")
               #else
            #       paste0("Z", outcome, "[", j_i, ", 1] * u[id", outcome, "[", j_i, "], ",
             #             RE_ind, "]"),
               paste0("Z", outcome, "[", j_i, ", 1] * u[id", outcome, "[", j_i, "], ",
                      RE_ind, "]", collapse = " + "),
               if (length(clm <- colmns_nHC)) {
                   #if (length(clm) > 1)
                #       paste0("\n", myt(4), " + dot_product(X", outcome, "[", j_i, ", colmns_nHC",
                  #           outcome,"], betas", outcome, "[colmns_nHC", outcome, "]);\n")
                 #  else
                   
                   paste0("\n", myt(4), " + ", 
                          paste0("X", outcome, "[", j_i, ", ", clm, "] * betas", outcome, 
                              "[", clm, "]", collapse = " + "), ";\n")
                   
               } else ";\n",
               myt(), "}\n")
    }
    paste0("\ntransformed parameters {\n", paste0(mapply(def_etas, outcomes), collapse = ""), 
           myt(), "matrix[n, n_RE] mu_u;\n",
           myt(), "for (i in 1:n) {\n",
           paste0(mapply(HC_part, rep(outcomes, sapply(colmns_HC, length)), 
                         colmns_HC2, nams, ncols, seq_along(colmns_HC2)), collapse = ""),
           myt(), "}\n",
           paste0(mapply(linpred_part, outcomes, colmns_HC, colmns_nHC, RE_inds), collapse = ""),
           "}\n")
}

##########################################################################################

model <- function (families, n_RE) {
    outcomes <- seq_along(families)
    RE_part <- paste0("\nmodel {\n",
                      if (n_RE > 1)
                          paste0(myt(),  "matrix[n_RE, n_RE] L_D;\n",
                                 myt(),  "L_D = diag_pre_multiply(L_var_D, L_corr_D);\n",
                                 myt(),  "L_var_D ~ student_t(3, 0, scale_diag_D);\n",
                                 myt(),  "L_corr_D ~ lkj_corr_cholesky(lkj_shape);\n",
                                 myt(),  "for (i in 1:n) {\n",
                                 myt(2), "u[i, ] ~ multi_normal_cholesky(mu_u[i, ], L_D);\n")
                      else
                          paste0(myt(), "D ~ student_t(3, 0, scale_diag_D);\n",
                                 myt(), "for (i in 1:n) {\n",
                                 myt(2), "u[i, ] ~ normal(mu_u[i, ], D);\n"),
                      myt(), "}\n")
    priors_part <- function (outcome, family) {
        paste0(myt(), "for (k", outcome, " in 1:ncx", outcome, ") {\n",
               myt(2), "betas", outcome, "[k", outcome, "] ~ normal(0.0, scale_betas", 
               outcome, ");\n",
               myt(), "}\n",
               if (family$family == "gaussian") paste0(myt(), "sigma", outcome, 
                                                       " ~ student_t(3, 0, scale_sigmas);\n"))
    }
    dist_part <- function (outcome, family) {
        if (family$family == "gaussian") {
            paste0(myt(), "y", outcome, " ~ normal(eta", outcome, ", sigma", outcome, ");\n")
        } else if (family$family == "binomial") {
            switch (family$link,
                "logit" = paste0(myt(), "y", outcome, " ~ bernoulli_logit(eta", outcome, ");\n"),
                "probit" = paste0(myt(), "y", outcome, " ~ bernoulli(Phi_approx(eta", outcome, "));\n"),
                "cloglog" = paste0(myt(), "y", outcome, " ~ bernoulli(inv_cloglog(eta", outcome, "));\n")
            )
        } else if (family$family == "poisson") {
            paste0(myt(), "y", outcome, " ~ poisson_log(eta", outcome, ");\n")
        }
    }
    paste0(RE_part, paste0(mapply(priors_part, outcomes, families), collapse = ""),
           paste0(mapply(dist_part, outcomes, families), collapse = ""), "}\n")
}

##########################################################################################

generated_quantities <- function (n_RE) {
    paste0("\ngenerated quantities {\n",
           if (n_RE > 1) paste0(myt(), "matrix[n_RE, n_RE] D;\n"),
           myt(), "matrix[n, n_RE] b;\n",
           if (n_RE > 1) paste0(myt(), "D = diag_pre_multiply(L_var_D, L_corr_D) * ",
                                "diag_pre_multiply(L_var_D, L_corr_D)';\n"),
           myt(), "b = u - mu_u;\n",
           "}\n")
}

##########################################################################################


#cat(data_part(families, lapply(colmns_HC, length), lapply(colmns_nHC, length), Data$n_RE),
#    parameters(families, Data$n_RE),
#    transformed_parameters(families, colmns_HC, colmns_nHC, RE_inds),
#    model(families, Data$n_RE),
#    generated_quantities(Data$n_RE))



