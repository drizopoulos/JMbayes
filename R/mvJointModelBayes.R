mvJointModelBayes <- function (mvglmerObject, survObject, timeVar,
                               Formulas = list(NULL), Interactions = list(NULL),
                               transFuns = NULL, priors = NULL, control = NULL, ...) {
    cl <- match.call()
    # control values
    con <- list(temps = 1.0, n_iter = 300, n_burnin = 1000,
                n_block = 50, n_thin = 300, target_acc = 0.234, c0 = 1, c1 = 0.8,
                eps1 = 1e-06, eps2 = 1e-05, eps3 = 1e04, adaptCov = FALSE,
                knots = NULL, ObsTimes.knots = TRUE,
                lng.in.kn = 15L, ordSpline = 4L, diff = 2L, speed_factor = 0.6,
                GQsurv = "GaussKronrod", GQsurv.k = 15L, seed = 1L,
                n_cores = max(1, parallel::detectCores() - 1), update_RE = TRUE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    # build survival data
    if (is.null(survObject$model)) {
        stop("'survObject' must be a 'coxph' or 'survreg' object fitted with argument 'model'",
             " set to TRUE.\n")
    }
    dataS <- survObject$model
    Terms <- attr(dataS, "terms")
    SurvInf <- model.response(dataS)
    typeSurvInf <- attr(SurvInf, "type")
    TimeVar <- all.vars(Terms)[1L]
    if (typeSurvInf == "right") {
        if (class(survObject)[1L] == 'survreg') {
            stop("Please refit the survival submodel using coxph().\n")
        }
        Time <- SurvInf[, "time"]
        Time[Time < 1e-04] <- 1e-04
        nT <- length(Time)
        event <- SurvInf[, "status"]
        LongFormat <- FALSE
        TimeLl <- rep(0.0, length(Time))
    }
    if (typeSurvInf == "counting") {
        if (class(survObject) == 'survreg') {
            stop("Please refit the survival submodel using coxph().\n")
        }
        if (is.null(survObject$model$cluster)) {
            stop("you need to refit the Cox and include in the right hand side of the ",
                 "formula the 'cluster()' function using as its argument the subjects' ",
                 "id indicator. These ids need to be the same as the ones used to fit ",
                 "the mixed effects model.\n")
        }
        idT <- survObject$model$cluster
        LongFormat <- length(idT) > length(unique(idT))
        TimeL <- TimeLl <- SurvInf[, "start"]
        fidT <- factor(idT, levels = unique(idT))
        TimeL <- tapply(TimeL, fidT, head, n = 1)
        anyLeftTrunc <- any(TimeL > 1e-07)
        TimeR <- SurvInf[, "stop"]
        TimeR[TimeR < 1e-04] <- 1e-04
        Time <- tapply(TimeR, fidT, tail, n = 1)
        nT <- length(Time)
        eventLong <- SurvInf[, "status"]
        event  <- c(tapply(eventLong, fidT, tail, n = 1))
        Terms <- drop.terms(Terms, attr(Terms,"specials")$cluster - 1,
                            keep.response = TRUE)
    }
    if (typeSurvInf == "interval") {
        Time1 <- SurvInf[, "time1"]
        Time2 <- SurvInf[, "time2"]
        Time <- Time1
        Time[Time2 != 1] <- Time2[Time2 != 1]
        TimeL <- Time1
        TimeL[Time2 == 1] <- 0.0
        Time[Time < 1e-04] <- 1e-04
        nT <- length(Time)
        event <- SurvInf[, "status"]
        LongFormat <- FALSE
        TimeLl <- rep(0.0, length(Time))
    }
    # Gauss-Kronrod components
    GQsurv <- if (con$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(con$GQsurv.k)
    wk <- GQsurv$wk
    sk <- GQsurv$sk
    K <- length(sk)
    P <- if (typeSurvInf == "counting") (Time - TimeL) / 2 else Time / 2
    st <- if (typeSurvInf == "counting") {
        outer(P, sk) + c(Time + TimeL) / 2
    } else {
        outer(P, sk + 1)
    }
    if (typeSurvInf == "interval") {
        P_int <- TimeL / 2
        st_int <- outer(P_int, sk + 1)
    }
    idGK <- rep(seq_len(nT), each = K)
    idGK_fast <- c(idGK[-length(idGK)] != idGK[-1L], TRUE)
    # knots baseline hazard
    kn <- if (is.null(con$knots)) {
        tt <- if (con$ObsTimes.knots) Time else Time[event == 1]
        pp <- quantile(tt, c(0.05, 0.95), names = FALSE)
        tail(head(seq(pp[1L], pp[2L], length.out = con$lng.in.kn), -1L), -1L)
    } else {
        con$knots
    }
    kn <- kn[kn < max(Time)]
    rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
    con$knots <- rr
    # build desing matrices for longitudinal process
    dataL <- mvglmerObject$data
    components <- mvglmerObject$components
    families <- mvglmerObject$families
    n_outcomes <- length(families)
    seq_n_outcomes <- seq_len(n_outcomes)
    idL <- components[paste0("id", seq_n_outcomes)]
    y <- components[paste0("y", seq_n_outcomes)]
    X <- components[paste0("X", seq_n_outcomes)]
    Z <- components[paste0("Z", seq_n_outcomes)]
    # constuct data set with the last repeated measurement per subject
    # evaluated at Time and at TimeL when interval censoring
    idVar <- components$idVar1
    if (is.null(dataL[[timeVar]])) {
        stop("variable '", timeVar, "' not in the data.frame extracted from 'mvglmerObject'.\n")
    }
    dataL <- dataL[order(dataL[[idVar]], dataL[[timeVar]]), ]
    last_rows <- function (data, ids) {
        fidVar <- factor(ids, levels = unique(ids))
        data[tapply(row.names(data), fidVar, tail, n = 1L), ]
    }
    dataL.id <- last_rows(dataL, dataL[[idVar]])
    dataL.id[[timeVar]] <- Time
    if (typeSurvInf == "interval") {
        dataL_int.id <- dataL.id
        dataL_int.id[[timeVar]] <- TimeL
    }
    # create the data set used for the calculation of the cumulative hazard
    # for the specified Gaussian quadrature points use the rows from the original
    # longitudinal data set that correspond to these points; this is to account for
    # time-varying covariates in the longitudinal submodel
    right_rows <- function (data, times, ids, Q_points) {
        fids <- factor(ids, levels = unique(ids))
        ind <- mapply(findInterval, split(Q_points, row(Q_points)), split(times, fids))
        ind[ind < 1] <- 1
        rownams_id <- split(row.names(data), fids)
        ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
        data[c(ind), ]
    }
    dataL.id2 <- right_rows(dataL, dataL[[timeVar]], dataL[[idVar]], st)
    dataL.id2[[timeVar]] <- c(t(st))
    if (typeSurvInf == "interval") {
        dataL_int.id2 <- right_rows(dataL, dataL[[timeVar]], dataL[[idVar]], st_int)
        dataL_int.id2[[timeVar]] <- c(t(st_int))
    }
    # create the survival data used for the calculation of the cumulative hazard
    # also create a merged data set from the longitudinal & survival submodels. The
    # latter is used to calculate interaction terms
    if (typeSurvInf == "right" || typeSurvInf == "interval") {
        idT <- dataS[[idVar]] <- unique(dataL[[idVar]])
    } else {
        if (!all(idT %in% unique(dataL[[idVar]]))) {
            stop("it seems there are some ids in the survival data set that cannot be ",
                 "found in the longitudinal data set.\n")
        }
        dataS[[idVar]] <- idT
    }
    dataS.id <- last_rows(dataS, dataS[[idVar]])
    dataS.id2 <- right_rows(dataS, TimeLl, idT, st)
    survVars_notin_long <- survVars_notin_long2 <- !names(dataS) %in% names(dataL)
    survVars_notin_long[names(dataS) == idVar] <- TRUE
    dataLS <- merge(dataL, dataS.id[survVars_notin_long], all = TRUE, sort = FALSE)
    dataLS.id <- merge(dataL.id, dataS.id[survVars_notin_long], by = idVar,
                       all = TRUE, sort = FALSE)
    dataLS.id[[TimeVar]] <- dataLS.id[[timeVar]]
    dataS.id2[["id2merge"]] <- paste(dataS.id2[[idVar]], round(c(t(st)), 8), sep = ":")
    dataL.id2[["id2merge"]] <- paste(dataL.id2[[idVar]], round(c(t(st)), 8), sep = ":")
    dataLS.id2 <- merge(dataL.id2, dataS.id2[survVars_notin_long2], by = "id2merge",
                        sort = FALSE, all = FALSE)
    dataLS.id2[[TimeVar]] <- dataLS.id2[[timeVar]]
    if (typeSurvInf == "interval") {
        dataS_int.id <- last_rows(dataS, dataS[[idVar]])
        dataS_int.id2 <- right_rows(dataS, TimeLl, idT, st_int)
        dataLS_int <- merge(dataL, dataS_int.id[survVars_notin_long], all = TRUE, 
                            sort = FALSE)
        dataLS_int.id <- merge(dataL_int.id, dataS_int.id[survVars_notin_long], by = idVar,
                           all = TRUE, sort = FALSE)
        uu <- runif(length(st_int))
        dataS_int.id2[["id2merge"]] <- paste(dataS_int.id2[[idVar]], uu, sep = ":")
        dataL_int.id2[["id2merge"]] <- paste(dataL.id2[[idVar]], uu, sep = ":")
        dataLS_int.id2 <- merge(dataL_int.id2, dataS_int.id2[survVars_notin_long2], 
                                by = "id2merge", sort = FALSE, all = FALSE)
    }
    # design matrices for the survival submodel, W1 is for the baseline hazard,
    # W2 for the baseline and external time-varying covariates
    W1 <- splines::splineDesign(con$knots, Time, ord = con$ordSpline)
    W1s <- splines::splineDesign(con$knots, c(t(st)), ord = con$ordSpline)
    W2 <- model.matrix(Terms, data = dataS.id)[, -1, drop = FALSE]
    W2s <- model.matrix(Terms, data = dataS.id2)[, -1, drop = FALSE]
    if (typeSurvInf == "interval") {
        W1_int <- splines::splineDesign(con$knots, TimeL, ord = con$ordSpline, 
                                        outer.ok = TRUE)
        W1s_int <- splines::splineDesign(con$knots, c(t(st_int)), ord = con$ordSpline, 
                                         outer.ok = TRUE)
        W2_int <- model.matrix(Terms, data = dataS_int.id)[, -1, drop = FALSE]
        W2s_int <- model.matrix(Terms, data = dataS_int.id2)[, -1, drop = FALSE]
    }
    extract_component <- function (component, fixed = TRUE) {
        components[grep(component, names(components), fixed = fixed)]
    }
    # expand the Formulas argument
    respVars <- unlist(extract_component("respVar"), use.names = FALSE)
    if (any(!names(Formulas) %in% respVars)) {
        stop("unknown names in the list provided in the 'Formulas' argument; as names ",
             "of the elements of this list you need to use the response variables from ",
             "the multivariate mixed model.\n")
    }
    # for outcomes not specified in Formulas use the value parameterization
    not_specified <- !respVars %in% names(Formulas)
    Formulas_ns <- rep(list("value"), length = sum(not_specified))
    names(Formulas_ns) <- respVars[not_specified]
    Formulas <- c(Formulas, Formulas_ns)
    Formulas <- Formulas[order(match(names(Formulas), respVars))]
    Formulas <- Formulas[!sapply(Formulas, is.null)]
    # extract the terms for the X and Z matrices to be used in creating the
    # corresponding design matrices
    TermsX <- extract_component("TermsX")
    X <- extract_component("^X[1-9]", FALSE)
    TermsZ <- extract_component("TermsZ")
    Z <- extract_component("^Z[1-9]", FALSE)
    names(TermsX) <- names(X) <- names(TermsZ) <- names(Z) <- respVars
    #
    which_value <- sapply(Formulas, function (x) any(x == "value"))
    names_which_value <- names(which_value)[which_value]
    replace_value <- function (termsx, x, termsz, z) {
        list(fixed = formula(termsx), indFixed = seq_len(ncol(x)),
             random = formula(termsz), indRandom = seq_len(ncol(z)),
             name = "value")
    }
    Formulas[which_value] <- mapply(replace_value, TermsX[names_which_value],
                                    X[names_which_value], TermsZ[names_which_value],
                                    Z[names_which_value], SIMPLIFY = FALSE)
    names_alphas <- function (Form) {
        name_term <- ifelse(!sapply(Form, is.list), "value", "extra")
        ind_extra <- name_term == "extra"
        user_names <- lapply(Form[ind_extra], "[[", "name")
        ind_usernames <- !sapply(user_names, is.null)
        name_term[ind_extra][ind_usernames] <- unlist(user_names[ind_usernames],
                                                      use.names = FALSE)
        out <- paste0(names(Form), "_", name_term)
        which_dupl <- unique(out[duplicated(out)])
        replc <- function (x) {
            paste(x, seq_along(x), sep = ".")
        }
        replacement <- unlist(lapply(which_dupl, function (dbl) replc(out[out == dbl])),
                              use.names = FALSE)
        out[out %in% which_dupl] <- replacement
        out
    }
    outcome <- match(names(Formulas), respVars)
    names(Formulas) <- names_alphas(Formulas)
    build_model_matrix <- function (input_terms, dataOrig, data, which) {
        out <- vector("list", length(input_terms))
        for (i in seq_along(input_terms)) {
            tr <- terms(input_terms[[i]][[which]], data = dataOrig)
            out[[i]] <- model.matrix(tr, model.frame(tr, data = data, na.action = NULL))
        }
        out
    }
    XX <- build_model_matrix(Formulas, dataL, dataL.id, "fixed")
    XXs <- build_model_matrix(Formulas, dataL, dataL.id2, "fixed")
    ZZ <- build_model_matrix(Formulas, dataL, dataL.id, "random")
    ZZs <- build_model_matrix(Formulas, dataL, dataL.id2, "random")
    if (typeSurvInf == "interval") {
        XXs_int <- build_model_matrix(Formulas, dataL, dataL_int.id2, "fixed")
        ZZs_int <- build_model_matrix(Formulas, dataL, dataL_int.id2, "random")
    }
    possible_names <- unique(c(respVars, paste0(respVars, "_value"), names(Formulas)))
    if (any(!names(Interactions) %in% possible_names)) {
        stop("unknown names in the list provided in the 'Interactions' argument; as names ",
             "of the elements of this list you need to use the response variables from ",
             "the multivariate mixed model or the induced names from the 'Formulas' ",
             "argument; these are: ", paste(names(Formulas), collapse = ", "), ".\n")
    }
    # replace names from Interactions to match the ones from Formulas
    ind_nams <- unlist(lapply(respVars, function (nam) which(names(Interactions) == nam)))
    names(Interactions)[ind_nams] <- paste0(names(Interactions)[ind_nams], "_value")
    if (any(duplicated(names(Interactions)))) {
        stop("duplicated names in argument 'Interactions'; check the help page.\n")
    }
    not_specified <- !names(Formulas) %in% names(Interactions)
    Interactions_ns <- rep(list(~ 1), length = sum(not_specified))
    names(Interactions_ns) <- names(Formulas)[not_specified]
    Interactions <- c(Interactions, Interactions_ns)
    Interactions <- Interactions[order(match(names(Interactions), names(Formulas)))]
    Interactions <- Interactions[!sapply(Interactions, is.null)]
    U <- lapply(Interactions, function (form) {
        model.matrix(terms(form, data = dataLS), data = dataLS.id)
    })
    Us <- lapply(Interactions, function (form) {
        model.matrix(terms(form, data = dataLS), data = dataLS.id2)
    })
    if (typeSurvInf == "interval") {
        Us_int <- lapply(Interactions, function (form) {
            model.matrix(terms(form, data = dataLS_int), data = dataLS_int.id2)
        })
    }
    id <- lapply(seq_len(n_outcomes), function (i) seq_len(nT))
    ids <- rep(list(idGK), n_outcomes)
    # extract fixed and random effects
    betas <- mvglmerObject$mcmc[grep("betas", names(mvglmerObject$mcmc), fixed = TRUE)]
    colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
    RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                      sq = components[grep("ncz", names(components), fixed = TRUE)],
                      incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                      SIMPLIFY = FALSE)
    bb <- mvglmerObject$mcmc$b
    b <- lapply(RE_inds, function (ind) bb[, , ind, drop = FALSE])
    inv_D <- mvglmerObject$mcmc[grep("inv_D", names(mvglmerObject$mcmc), fixed = TRUE)]
    sigmas <- vector("list", n_outcomes)
    if (any(which_gaussian <- sapply(families, "[[", "family") == "gaussian")) {
        sigmas[which_gaussian] <- mvglmerObject$mcmc[grep("sigma", 
                                                          names(mvglmerObject$mcmc), 
                                                          fixed = TRUE)]
    }
    # create design matrix long for relative risk model
    trans_Funs <- rep("identity", length(XX))
    names(trans_Funs) <- names(Formulas)
    if (!is.null(transFuns)) {
        if (is.null(names(transFuns)) || !all(names(transFuns) %in% names(trans_Funs))) {
            stop("unknown names in 'transFuns'; valid names are the ones induced by the 'Formulas' ",
                 "argument; these are:\n", paste(names(Formulas), collapse = ", "))
        }
        valid_funs <- c("identity", "expit", "exp", "log", "log2", "log10", "sqrt")
        if (!is.character(transFuns) || !all(transFuns %in% valid_funs)) {
            stop("invalid functions names in 'transFuns'; the functions currently supported are: ", 
                 paste(valid_funs, collapse = ", "), ".\nThese should be provided as character vector.")
        }
        trans_Funs[names(transFuns)] <- transFuns
    }
    indFixed <- lapply(Formulas, "[[", "indFixed")
    indRandom <- lapply(Formulas, "[[", "indRandom")
    RE_inds2 <- mapply(function (ind, select) ind[select], RE_inds[outcome], indRandom,
                       SIMPLIFY = FALSE)
    Xbetas_calc <- function (X, betas, index = NULL, outcome) {
        n <- length(X)
        out <- vector("list", n)
        for (i in seq_len(n)) {
            out[[i]] <- if (is.null(index)) {
                c(X[[i]] %*% betas[[i]])
            } else {
                betas_i <- betas[[outcome[i]]]
                c(X[[i]] %*% betas_i[index[[i]]])
            }
        }
        out
    }
    get_fun <- function (f) {
        if (f == "identity") {
            function (x) x 
        } else if (f == "expit") {
            function (x) exp(x) / (1 + exp(x))
        } else get(f, mode = "function")
    }
    designMatLong <- function (X, betas, Z, b, id, outcome, indFixed, indRandom, U,
                               trans_Funs) {
        n <- length(X)
        cols <- sapply(U, ncol)
        cols_inds <- cbind(c(1, head(cumsum(cols) + 1, -1)), cumsum(cols))
        n_out <- sum(cols)
        col_inds_out <- vector("list", n)
        out <- matrix(0, nrow(X[[1]]), n_out)
        for (i in seq_len(n)) {
            ii <- outcome[i]
            iii <- col_inds_out[[i]] <- seq(cols_inds[i, 1], cols_inds[i, 2])
            X_i <- X[[i]]
            betas_i <- betas[[ii]][indFixed[[i]]]
            Z_i <- Z[[i]]
            b_i <- as.matrix(b[[ii]])[id[[ii]], indRandom[[i]], drop = FALSE]
            Fun <- get_fun(trans_Funs[i])
            out[, iii] <- U[[i]] * Fun(c(X_i %*% betas_i) + rowSums(Z_i * b_i))
        }
        attr(out, "col_inds") <- col_inds_out
        out
    }
    postMean_betas <- lapply(betas, colMeans, na.rm = TRUE)
    postMean_b <- lapply(b, function (m) apply(m, 2:3, mean, na.rm = TRUE))
    postMean_inv_D <- lapply(inv_D, function (m) apply(m, 2:3, mean, na.rm = TRUE))
    mean_null <- function (x) if (is.null(x)) as.numeric(NA) else mean(x)
    postMean_sigmas <- lapply(sigmas, mean_null)
    Xbetas <- Xbetas_calc(X, postMean_betas)
    XXbetas <- Xbetas_calc(XX, postMean_betas, indFixed, outcome)
    XXsbetas <- Xbetas_calc(XXs, postMean_betas, indFixed, outcome)
    if (typeSurvInf == "interval") {
        XXsbetas_int <- Xbetas_calc(XXs_int, postMean_betas, indFixed, outcome)
    }
    fams <- sapply(families, "[[", "family")
    links <- sapply(families, "[[", "link")
    idL2 <- lapply(idL, function (x) {
        x <- c(x[-length(x)] != x[-1L], TRUE)
        which(x) - 1
    })
    Wlong <- designMatLong(XX, postMean_betas, ZZ, postMean_b, id, outcome,
                           indFixed, indRandom, U, trans_Funs)
    Wlongs <- designMatLong(XXs, postMean_betas, ZZs, postMean_b, ids, outcome,
                            indFixed, indRandom, Us, trans_Funs)
    if (typeSurvInf == "interval") {
        Wlongs_int <- designMatLong(XXs_int, postMean_betas, ZZs_int, postMean_b, ids, 
                                    outcome, indFixed, indRandom, Us_int)
    }
    # priors
    DD <- diag(ncol(W1))
    Tau_Bs_gammas <- crossprod(diff(DD, differences = con$diff)) + 1e-06 * DD
    find_td_cols <- function (x) grep('tve(', colnames(x), fixed = TRUE)
    td_cols <- lapply(U, find_td_cols)
    Tau_alphas <- lapply(U, function (x) 0.01 * diag(NCOL(x)))
    pen_matrix <- function (td_cols, Tau_alphas) {
        if (ncol(Tau_alphas) == length(td_cols)) {
            DD <- diag(ncol(Tau_alphas))
            crossprod(diff(DD, differences = con$diff)) + 1e-06 * DD
        } else Tau_alphas
    }
    Tau_alphas <- bdiag(mapply(pen_matrix, td_cols, Tau_alphas, SIMPLIFY = FALSE))
    which_td <- as.logical(sapply(td_cols, length))
    td_cols <- mapply(seq, from = c(1, head(cumsum(sapply(U, ncol)), -1) + 1),
                      to = cumsum(sapply(U, ncol)), SIMPLIFY = FALSE)[which_td]
    prs <- list(mean_Bs_gammas = rep(0, ncol(W1)), Tau_Bs_gammas = Tau_Bs_gammas,
                mean_gammas = rep(0, ncol(W2)), Tau_gammas = 0.01 * diag(ncol(W2)),
                mean_alphas = rep(0, ncol(Wlong)), Tau_alphas = Tau_alphas,
                td_cols = unname(td_cols),
                rank_Tau_td_alphas = if (length(td_cols)) length(td_cols[[1]]) else 0,
                A_tau_Bs_gammas = 1, B_tau_Bs_gammas = 0.01, 
                rank_Tau_Bs_gammas = qr(Tau_Bs_gammas)$rank,
                A_phi_Bs_gammas = 1, B_phi_Bs_gammas = 0.01, shrink_Bs_gammas = FALSE,
                A_tau_gammas = 0.1, B_tau_gammas = 0.1, rank_Tau_gammas = ncol(W2),
                A_phi_gammas = 1, B_phi_gammas = 0.01, shrink_gammas = FALSE,
                A_tau_alphas = 0.1, B_tau_alphas = 0.1, rank_Tau_alphas = ncol(Wlong),
                A_phi_alphas = 1, B_phi_alphas = 0.01, shrink_alphas = FALSE)
    if (!is.null(priors)) {
        lngths <- lapply(prs[(nam.prs <- names(priors))], length)
        if (!is.list(priors) || !isTRUE(all.equal(lngths, lapply(priors, length)))) {
            warning("'priors' is not a list with elements numeric vectors of appropriate ",
                    "length; default priors are used instead.\n")
        } else {
            prs[nam.prs] <- priors
        }
    }
    if (mvglmerObject$engine == "JAGS") {
        tau_betas <- mvglmerObject$priors[grep("tau_betas", names(mvglmerObject$priors), 
                                               fixed = TRUE)]
        prs$Tau_betas <- diag(rep(unlist(tau_betas, use.names = FALSE), 
                                  sapply(betas, ncol)))
    } else {
        scale_betas <- mvglmerObject$priors[grep("scale_betas", names(mvglmerObject$priors), 
                                               fixed = TRUE)]
        prs$Tau_betas <- diag(rep(1 / unlist(scale_betas, use.names = FALSE)^2, 
                                  sapply(betas, ncol)))
    }
    prs$priorK_D <- mvglmerObject$priors$priorK_D
    # Data passed to the MCMC
    Data <- list(y = y, Xbetas = Xbetas, X = X, Z = Z, RE_inds = RE_inds,
                 RE_inds2 = RE_inds2, idL = idL, idL2 = idL2, sigmas = postMean_sigmas,
                 invD = postMean_inv_D[[1]], fams = fams, links = links, Time = Time,
                 event = event, idGK_fast = which(idGK_fast) - 1, W1 = W1, W1s = W1s,
                 event_colSumsW1 = colSums(event * W1), W2 = W2, W2s = W2s,
                 event_colSumsW2 = if (ncol(W2)) colSums(event * W2),
                 Wlong = Wlong, Wlongs = Wlongs,
                 event_colSumsWlong = colSums(event * Wlong),
                 U = U, Us = Us, col_inds = attr(Wlong, "col_inds"),
                 row_inds_U = seq_len(nrow(Wlong)), row_inds_Us = seq_len(nrow(Wlongs)),
                 XXbetas = XXbetas, XXsbetas = XXsbetas, XX = XX, XXs = XXs, ZZ = ZZ,
                 ZZs = ZZs, P = P[ids[[1]]], w = rep(wk, nT),
                 Pw = P[ids[[1]]] * rep(wk, nT), idT = id[outcome], idTs = ids[outcome],
                 outcome = outcome, indFixed = indFixed, indRandom = indRandom,
                 trans_Funs = trans_Funs)
    if (typeSurvInf == "interval") {
        Data <- c(Data, list(W1s_int = W1s_int, W2s_int = W2s_int, Wlongs_int = Wlongs_int,
                             Us_int = Us_int, XXsbetas_int = XXsbetas_int,
                             XXs_int = XXs_int, ZZs_int = ZZs_int,
                             P_int = P_int[ids[[1]]], 
                             Pw_int = P_int[ids[[1]]] * rep(wk, nT)))
    }
    # initial values
    inits <- list(Bs_gammas = rep(0, ncol(W1)), tau_Bs_gammas = 200, phi_Bs_gammas = rep(1, ncol(W1)),
                  gammas = rep(0, ncol(W2)), tau_gammas = 1, phi_gammas = rep(1, ncol(W2)),
                  alphas = rep(0, ncol(Wlong)), tau_alphas = 1, phi_alphas = rep(1, ncol(Wlong)),
                  tau_td_alphas = rep(200, length(td_cols)))
    inits2 <- marglogLik2(inits[c("Bs_gammas", "gammas", "alphas", "tau_Bs_gammas")],
                          Data, prs, fixed_tau_Bs_gammas = TRUE)
    inits[names(attr(inits2, "inits"))] <- attr(inits2, "inits")
    Cvs <- attr(inits2, "Covs")
    nRE <- sum(sapply(Z, ncol))
    Cvs$b <- array(0.0, c(nRE, nRE, nT))
    for (i in seq_len(nT)) Cvs$b[, , i] <- chol(var(bb[, i, ]))
    if (is.null(Cvs$gammas)) Cvs$gammas <- matrix(nrow = 0, ncol = 0)
    inits$b <- do.call("cbind", postMean_b)
    scales <- list(b = rep(5.76 / nRE, nT), Bs_gammas = 5.76 / ncol(W1),
                   gammas = 5.76 / ncol(W2), alphas = 5.76 / ncol(Wlong))
    sampl <- function (x, m) {
        lapply(x, function (obj) {
            d <- dim(obj)
            if (is.null(d)) {
                obj[m]
            } else {
                if (is.matrix(obj)) obj[m, ] else obj[m, , ]
            }

        })
    }
    runParallel <- function (block, betas, b, sigmas, inv_D, inits, data, priors,
                             scales, Covs, control, interval_cens) {
        M <- length(block)
        LogLiks <- numeric(M)
        out <- vector("list", M)
        new_scales <- vector("list", M)
        inits_Laplace <- inits[c("Bs_gammas", "gammas", "alphas", "tau_Bs_gammas")]
        inits_Laplace[["tau_Bs_gammas"]] <- log(inits_Laplace[["tau_Bs_gammas"]])
        inits_Laplace[["b"]] <- NULL
        any_gammas <- as.logical(length(priors[["mean_gammas"]]))
        set.seed(control$seed)
        on.exit(rm(list = ".Random.seed", envir = globalenv()))
        for (i in seq_len(M)) {
            ii <- block[i]
            if (control$update_RE) {
                betas. <- sampl(betas, ii)
                data$Xbetas <- Xbetas_calc(data$X, betas.)
                outcome <- data$outcome
                indFixed <- data$indFixed
                data$XXbetas <- Xbetas_calc(data$XX, betas., indFixed, outcome)
                data$XXsbetas <- Xbetas_calc(data$XXs, betas., indFixed, outcome)
                if (interval_cens) {
                    data$XXsbetas_int <- Xbetas_calc(data$XXs_int, betas., indFixed, outcome)
                }
                data$sigmas <- sampl(sigmas, ii)
                data$invD <- as.matrix(sampl(inv_D, ii)[[1]])
                oo <- lap_rwm_C(inits, data, priors, scales, Covs, control, interval_cens)
                current_betas <- unlist(betas., use.names = FALSE)
                n_betas <- length(current_betas)
                pr_betas <- c(dmvnorm2(rbind(current_betas), rep(0, n_betas),
                                      priors$Tau_betas, logd = TRUE))
                pr_invD <- dwish(data$invD, diag(nrow(data$invD)),
                                 priors$priorK_D, log = TRUE)
                LogLiks[i] <- c(oo$logWeights) - pr_betas - pr_invD
                out[[i]] <- oo$mcmc
                new_scales[[i]] <- oo$scales$sigma
            } else {
                betas. <- sampl(betas, ii)
                b. <- sampl(b, ii)
                outcome <- data$outcome
                indFixed <- data$indFixed
                indRandom <- data$indRandom
                data$Wlong <- designMatLong(data$XX, betas., data$ZZ, b., data$idT,
                                            outcome, indFixed, indRandom, data$U)
                data$Wlongs <- designMatLong(data$XXs, betas., data$ZZs, b., data$idTs,
                                             outcome, indFixed, indRandom, data$Us)
                data$event_colSumsWlong <- colSums(data$event * data$Wlong)
                LogLiks[i] <- marglogLik2(inits_Laplace, data, priors)
                oo <- if (any_gammas) {
                    lap_rwm_C_woRE(inits, data, priors, scales, Covs, control)
                } else {
                    lap_rwm_C_woRE_nogammas(inits, data, priors, scales, Covs, control)
                }
                out[[i]] <- oo$mcmc
                new_scales[[i]] <- oo$scales$sigma
            }
        }
        out <- lapply(unlist(out, recursive = FALSE), drop)
        new_scales <- lapply(unlist(new_scales, recursive = FALSE), drop)
        nams <- names(out)
        if (!is.null(out$b)) {
            b_out <- array(0.0, c(dim(as.matrix(out$b)), M))
            for (i in seq_len(M)) b_out[, , i] <- out[nams == "b"][[i]]
        } else b_out <- NULL
        out <- list("b" = b_out,
                    "Bs_gammas" = do.call("rbind", out[nams == "Bs_gammas"]),
                    "gammas" = if (any_gammas) do.call("rbind", out[nams == "gammas"]),
                    "alphas" = do.call("rbind", out[nams == "alphas"]),
                    "tau_Bs_gammas" = do.call("rbind", out[nams == "tau_Bs_gammas"]),
                    "tau_gammas" = if (any_gammas)do.call("rbind", out[nams == "tau_gammas"]),
                    "tau_alphas" = do.call("rbind", out[nams == "tau_alphas"]),
                    "tau_td_alphas" = do.call("rbind", out[nams == "tau_td_alphas"]),
                    "phi_Bs_gammas" = do.call("rbind", out[nams == "phi_Bs_gammas"]),
                    "phi_gammas" = if (any_gammas)do.call("rbind", out[nams == "phi_gammas"]),
                    "phi_alphas" = do.call("rbind", out[nams == "phi_alphas"]))
        out$LogLiks <- LogLiks
        nams <- names(scales)
        out$scales <- list("b" = do.call("rbind", new_scales[nams == "b"]),
                           "Bs_gammas" = do.call("rbind", new_scales[nams == "Bs_gammas"]),
                           "gammas" = if (any_gammas) do.call("rbind", new_scales[nams == "gammas"]),
                           "alphas" = do.call("rbind", new_scales[nams == "alphas"]) )
        out <- out[!sapply(out, is.null)]
        list(mcmc = out)
    }
    any_gammas <- ncol(W2)
    combine <- function(lis) {
        f <- function (lis, nam) {
            if (nam == "LogLiks") {
                lis <- unlist(lis, recursive = FALSE)
                nam <- paste0("mcmc.", nam)
                unname(do.call("c", lis[names(lis) == nam]))
            } else {
                lis <- unlist(lis[names(lis) == "mcmc"], recursive = FALSE)
                nam <- paste0("mcmc.", nam)
                if (nam == "mcmc.b")
                    abind(lis[names(lis) == nam])
                else
                    unname(do.call("rbind", lis[names(lis) == nam]))
            }
        }
        list("b" = f(lis, "b"),
             "Bs_gammas" = f(lis, "Bs_gammas"), "tau_Bs_gammas" = f(lis, "tau_Bs_gammas"),
             "phi_Bs_gammas" = f(lis, "phi_Bs_gammas"),
             "gammas" = if (any_gammas) f(lis, "gammas"),
             "tau_gammas" = if (any_gammas) f(lis, "tau_gammas"),
             "phi_gammas" = if (any_gammas) f(lis, "phi_gammas"),
             "alphas" = f(lis, "alphas"), "tau_alphas" = f(lis, "tau_alphas"),
             "tau_td_alphas" = f(lis, "tau_td_alphas"),
             "phi_alphas" = f(lis, "phi_alphas"),
             "LogLiks" = f(lis, "LogLiks"))
    }
    # We first split the number of iterations in blocks according to the number of
    # processors; the first block is split again according to the number of processors;
    # first we run in parallel the first block; we the update the scales, and following
    # we run the rest of the original blocks
    M <- nrow(betas[[1L]])
    blocks <- split(seq_len(M),
                    rep(seq_len(con$n_cores + 1L), each = ceiling(M / (con$n_cores + 1L)),
                        length.out = M))
    block1 <- split(blocks[[1L]],
                    rep(seq_len(con$n_cores), each = ceiling(length(blocks[[1L]]) / con$n_cores),
                        length.out = length(blocks[[1L]])))
    blocks <- blocks[-1L]
    elapsed_time <- system.time({
        cluster <- makeCluster(con$n_cores)
        registerDoParallel(cluster)
        out1 <- foreach(i = block1, .packages = "JMbayes", .combine = c) %dopar% {
            runParallel(i, betas, b, sigmas, inv_D, inits, Data, prs, scales, Cvs, con,
                        typeSurvInf == "interval")
        }
        stopCluster(cluster)
        if (con$speed_factor < 1) {
            calc_new_scales <- function (parm) {
                if (parm == "b") {
                    apply(do.call("rbind", lapply(new_scales, "[[", "b")), 2L, median)
                } else {
                    median(do.call("c", lapply(new_scales, "[[", parm)))
                }
            }
            new_scales <- lapply(out1, "[[", "scales")
            new_scales <- list("b" = if (con$update_RE) calc_new_scales("b"),
                               "Bs_gammas" = calc_new_scales("Bs_gammas"),
                               "gammas" = if (any_gammas) calc_new_scales("gammas") else Inf,
                               "alphas" = calc_new_scales("alphas"))
            con$n_burnin <- round(ceiling(con$speed_factor * 
                                              con$n_burnin / con$n_block) * con$n_block)
        } else {
            new_scales <- scales
        }
        cluster <- makeCluster(con$n_cores)
        registerDoParallel(cluster)
        out <- foreach(i = blocks, .packages = "JMbayes", .combine = c) %dopar% {
            runParallel(i, betas, b, sigmas, inv_D, inits, Data, prs, new_scales, Cvs, 
                        con, typeSurvInf == "interval")
        }
        stopCluster(cluster)
        out <- c(out1, out)
    })["elapsed"]
    # collect the results
    mcmc <- mvglmerObject$mcmc
    keep <- unlist(sapply(c("betas", "sigma", "D"), grep, x = names(mcmc), fixed = TRUE))
    mcmc <- c(mcmc[keep], combine(out))
    # set the appropriate names
    colnames(mcmc$Bs_gammas) <- paste0("bs", seq_len(ncol(mcmc$Bs_gammas)))
    colnames(mcmc$gammas) <- colnames(W2)
    get_U_colnames <- unlist(lapply(U, function (u)
        gsub("(Intercept)", "", colnames(u), fixed = TRUE)))
    nam_alph <- names(U)
    trans_Funs <- trans_Funs[names(U)]
    ind <- trans_Funs != "identity"
    nam_alph[ind] <- paste0(trans_Funs[ind], "(", nam_alph, ")")
    colnames(mcmc$alphas) <- paste0(rep(nam_alph, sapply(U, ncol)),
                                    ifelse(get_U_colnames == "", "", ":"),
                                    get_U_colnames)
    # caclulate summaries from the Monte Carlo samples
    summary_fun <- function (FUN, ...) {
        fun <- function (x, ...) {
            res <- try(FUN(x, ...), silent = TRUE)
            if (!inherits(res, "try-error"))  res else NA
        }
        out <- lapply(mcmc, function (x) {
            d <- dim(x)
            if (!is.null(d) && length(d) > 1) {
                dd <- if (length(d) == 2) 2L else if (d[1L] > d[3L]) c(2L, 3L)
                else c(1L, 2L)
                apply(x, dd, fun, ...)
            } else if (!is.null(x)) {
                fun(x, ...)
            }
        })
        out[!sapply(out, is.null)]
    }
    stand <- function (x) {
        n <- length(x)
        upp <- max(x, na.rm = TRUE) + log(n)
        w <- exp(x - upp)
        w / sum(w)
    }
    LogLiks <- combine(out)$LogLiks
    weights <- stand(LogLiks)
    wmean <- function (x, weights, na.rm = FALSE) sum(x * weights, na.rm = na.rm)
    wsd <- function (x, weights) sqrt(drop(cov.wt(cbind(x[!is.na(x)]), wt = weights)$cov))
    # Extract the results
    res <- list(call = cl, mcmc = mcmc, weights = weights,
                mcmc_info = list(
                    elapsed_mins = elapsed_time / 60, 
                    n_burnin = con$n_burnin, 
                    n_iter = con$n_iter + con$n_burnin, n_thin = con$n_thin,
                    priors = prs
                ),
                statistics = list(
                    postMeans = summary_fun(mean, na.rm = TRUE),
                    postwMeans = summary_fun(wmean, weights = weights, na.rm = TRUE),
                    postModes = summary_fun(modes),
                    EffectiveSize = summary_fun(effectiveSize),
                    StDev = summary_fun(sd, na.rm = TRUE),
                    wStDev = summary_fun(wsd, weights = weights),
                    StErr = summary_fun(stdErr),
                    CIs = summary_fun(quantile, probs = c(0.025, 0.975)),
                    wCIS = summary_fun(Hmisc::wtd.quantile, weights = weights, 
                                       probs = c(0.025, 0.975), type = "i/(n+1)", 
                                       na.rm = TRUE),
                    Pvalues = summary_fun(computeP)
                ),
                model_info = list(
                    families = families,
                    timeVar = timeVar,
                    Formulas = Formulas,
                    Interactions = Interactions,
                    RE_inds = RE_inds,
                    RE_inds2 = RE_inds2,
                    transFuns = trans_Funs,
                    mvglmer_components = components,
                    coxph_components = list(data = dataS, Terms = Terms, Time = Time, 
                                            event = event),
                    functions = list(build_model_matrix = build_model_matrix,
                                     last_rows = last_rows, right_rows = right_rows,
                                     Xbetas_calc = Xbetas_calc, 
                                     designMatLong = designMatLong,
                                     get_fun = get_fun)
                ),
                control = con)
    class(res) <- "mvJMbayes"
    res
}


