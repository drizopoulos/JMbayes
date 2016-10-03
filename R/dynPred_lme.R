extract_lmeComponents <- function (lmeObject, timeVar) {
    data <- lmeObject$data
    formYx <- formula(lmeObject)
    mfX <- model.frame(terms(formYx), data = data)
    TermsX <- attr(mfX, "terms")
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    idVar <- names(lmeObject$modelStruct$reStruct)
    betas <- fixef(lmeObject)
    sigma <- lmeObject$sigma
    D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", sigma^2)[[1]]
    V <- vcov(lmeObject)
    times_orig <- data[[timeVar]]
    times_orig <- times_orig[!is.na(times_orig)]
    out <- list(formYx = formYx, TermsX = TermsX, formYz = formYz, TermsZ = TermsZ, 
                idVar = idVar, betas = betas, sigma = sigma, D = D, V = V,
                times_orig = times_orig)
    class(out) <- "lmeComponents"
    out
}

IndvPred_lme <- function (lmeObject, newdata, timeVar, times = NULL, M = 200L,
                          interval = c("confidence", "prediction"),
                          level = 0.95, return_data = FALSE, seed = 1L) {
    if (!inherits(lmeObject, "lme") && !inherits(lmeObject, "lmeComponents"))
        stop("Use only with 'lme' or 'lmeComponents' objects.\n")
    interval <- match.arg(interval)
    if (inherits(lmeObject, "lme")) {
        data <- lmeObject$data
        formYx <- formula(lmeObject)
        mfX <- model.frame(terms(formYx), data = data)
        TermsX <- attr(mfX, "terms")
        formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
        mfZ <- model.frame(terms(formYz), data = data)
        TermsZ <- attr(mfZ, "terms")
        idVar <- names(lmeObject$modelStruct$reStruct)
        betas <- fixef(lmeObject)
        sigma <- lmeObject$sigma
        D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", sigma^2)[[1]]
        V <- vcov(lmeObject)
        times_orig <- data[[timeVar]]
        times_orig <- times_orig[!is.na(times_orig)]
    } else {
        formYx <- lmeObject$formYx
        TermsX <- lmeObject$TermsX
        formYz <- lmeObject$formYz
        TermsZ <- lmeObject$TermsZ 
        idVar <- lmeObject$idVar
        betas <- lmeObject$betas
        sigma <- lmeObject$sigma
        D <- lmeObject$D
        V <- lmeObject$V
        times_orig <- lmeObject$times_orig
    }
    mfX_new <- model.frame(TermsX, data = newdata)
    X_new <- model.matrix(formYx, mfX_new)
    mfZ_new <- model.frame(TermsZ, data = newdata)
    Z_new <- model.matrix(formYz, mfZ_new)
    y_new <- model.response(mfX_new, "numeric")
    if (length(idVar) > 1)
        stop("the current version of the function only works with a single grouping variable.\n")
    if (is.null(newdata[[idVar]]))
        stop("subject id variable not in newdata.")
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    n <- length(unique(id))
    ######################################################################################
    modes <- matrix(0.0, n, ncol(Z_new))
    post_vars <- DZtVinv <- vector("list", n)
    for (i in seq_len(n)) {
        id_i <- id == i
        X_new_id <- X_new[id_i, , drop = FALSE]
        Z_new_id <- Z_new[id_i, , drop = FALSE]
        Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + sigma^2 * diag(sum(id_i)))
        DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
        modes[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - X_new_id %*% betas))
        t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
        t2 <- DZtVinv[[i]] %*% X_new_id %*% V %*% 
            crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
        post_vars[[i]] <- D - t1 + t2
    }
    fitted_y <- c(X_new %*% betas) + rowSums(Z_new * modes[id, , drop = FALSE])
    ######################################################################################
    if (is.null(times) || !is.numeric(times)) {
        times <- seq(min(times_orig), max(times_orig), length.out = 100)
    }
    newdata_pred <- newdata[tapply(row.names(newdata), id, tail, n = 1), ]
    last_time <- newdata_pred[[timeVar]]
    times_to_pred <- lapply(last_time, function (t) times[times > t])
    id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
    newdata_pred <- newdata_pred[id_pred, ]
    newdata_pred[[timeVar]] <- unlist(times_to_pred)
    mfX_new_pred <- model.frame(TermsX, data = newdata_pred)
    X_new_pred <- model.matrix(formYx, mfX_new_pred)
    mfZ_new_pred <- model.frame(TermsZ, data = newdata_pred)
    Z_new_pred <- model.matrix(formYz, mfZ_new_pred)
    predicted_y <- c(X_new_pred %*% betas) + 
        rowSums(Z_new_pred * modes[id_pred, , drop = FALSE])
    set.seed(seed)
    betas_M <- MASS::mvrnorm(M, betas, V)
    modes_fun <- function (betas) {
        t(mapply("%*%", DZtVinv, split(y_new - X_new %*% betas, id)))
    }
    modes_M <- lapply(split(betas_M, row(betas_M)), modes_fun)
    matrix_row <- function (m, i) m[i, , drop = FALSE]
    modes_M <- lapply(seq_len(n), function (i) t(sapply(modes_M, matrix_row, i = i)))
    b_M <- modes_M
    for (i in seq_len(n)) {
        b_M[[i]] <- t(apply(modes_M[[i]], 1, mvrnorm, n = 1, Sigma = post_vars[[i]]))
    }
    n_pred <- length(predicted_y)
    sampled_y <- matrix(0.0, n_pred, M)
    for (m in seq_len(M)) {
        betas_m <- betas_M[m, ]
        b_m <- t(sapply(b_M, function (x) x[m, ]))
        mean_m <- c(X_new_pred %*% betas_m) + 
            rowSums(Z_new_pred * b_m[id_pred, , drop = FALSE])
        sampled_y[, m] <- if (interval == "confidence") mean_m 
            else rnorm(n_pred, mean_m, lmeObject$sigma)
    }
    low <- apply(sampled_y, 1, quantile, probs = (1 - level) / 2)
    upp <- apply(sampled_y, 1, quantile, probs = 1 - (1 - level) / 2)
    rm(list = ".Random.seed", envir = globalenv())
    if (!return_data) {
        list(times_to_pred = times_to_pred, predicted_y = predicted_y, 
             low = low, upp = upp)
    } else {
        out_data <- rbind(newdata, newdata_pred)
        out_data$pred <- c(fitted_y, predicted_y)
        out_data$low <- c(rep(NA, length(fitted_y)), low)
        out_data$upp <- c(rep(NA, length(fitted_y)), upp)
        out_data[order(out_data[[idVar]], out_data[[timeVar]]), ]
    }
}