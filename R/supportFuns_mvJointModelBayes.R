build_model_matrix <- function (input_terms, dataOrig, data, which) {
    out <- vector("list", length(input_terms))
    for (i in seq_along(input_terms)) {
        MF <- model.frame.default(terms(input_terms[[i]][[which]]), dataOrig)
        tr <- terms(MF)
        out[[i]] <- model.matrix(tr, model.frame(tr, data = data, na.action = NULL))
    }
    out
}

last_rows <- function (data, ids) {
    fidVar <- factor(ids, levels = unique(ids))
    data[tapply(row.names(data), fidVar, tail, n = 1L), ]
}

right_rows <- function (data, times, ids, Q_points) {
    fids <- factor(ids, levels = unique(ids))
    if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
    ind <- mapply(findInterval, Q_points, split(times, fids))
    ind[ind < 1] <- 1
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
}

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
