extractFrames <- function (formula, data) {
    Terms <- terms(formula)
    term_labels <- attr(Terms, "term.labels")
    which_RE <- grep("|", term_labels, fixed = TRUE)
    namesVars <- all.vars(formula)
    respVar <- as.character(formula)[2L]
    # Fixed Effects
    formYx <- paste(term_labels[-which_RE], collapse = " + ")
    formYx <- as.formula(paste(respVar, "~", formYx))
    TermsX <- terms(formYx, data = data)
    mfX <- model.frame(TermsX, data)
    X <- model.matrix(TermsX, data)
    # Random Effects
    spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
    idVar <- spl[2L]
    data <- data[complete.cases(data[namesVars]), ]
    id <- data[[idVar]]
    id <- match(id, unique(id))
    formYz <- paste(spl[1], collapse = " + ")
    formYz <- as.formula(paste(respVar, "~", formYz))
    TermsZ <- terms(formYz, data = data)
    Z <- model.matrix(TermsZ, data)
    # response variable
    y <- model.response(mfX)
    if (is.factor(y))
        y <- as.vector(unclass(y) - 1)
    # hierarchical centering
    find_positions <- function (nams1, nams2) {
        nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
        vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
                  glob2rx(paste0("*:", nams1)))
        out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
        out
    }
    check_td <- function (x, id) {
        !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
    }
    has_interceptX <- attr(TermsX, "intercept")
    has_interceptZ <- attr(TermsZ, "intercept")
    performHC <- has_interceptX && (has_interceptX == has_interceptZ)
    if (performHC) {
        terms.labs_X <- attr(TermsX, "term.labels")
        terms.labs_Z <- attr(TermsZ, "term.labels")
        # check for time-varying covariates
        timeTerms <- if (length(terms.labs_Z)) 
            grep(terms.labs_Z[1L], colnames(X), fixed = TRUE)
        which_td <- unname(which(apply(X, 2, check_td, id = id)))
        all_TDterms <- unique(c(timeTerms, which_td))
        baseline <- seq_len(ncol(X))[-all_TDterms]
        ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions,
                                               nams2 = colnames(X)))
        ind_colmns2 <- seq_len(ncol(X))
        ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
        data.id <- data[!duplicated(id), ]
        Xhc <- if (length(terms.labs_Z)) {
            #data.id[[terms.labs_Z[1L]]] <- 1
            mfHC <- model.frame(TermsX, data = data.id)
            which.timeVar <- grep(terms.labs_Z[1L], names(mfHC), fixed = TRUE)
            mfHC[which.timeVar] <- lapply(mfHC[which.timeVar],
                                          function (x) { x[] <- 1; x })
            model.matrix(formYx, mfHC)
        } else {
            model.matrix(formYx, model.frame(TermsX, data = data.id))
        }
    }
    # extract results
    list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
         id = id, y = y, X = X, Z = Z, TermsX = TermsX,
         TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
         Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
         ncx = ncol(X), ncz = ncol(Z))
}

abind <- function (..., along = N, rev.along = NULL, new.names = NULL,
          force.array = TRUE, make.names = use.anon.names, use.anon.names = FALSE,
          use.first.dimnames = FALSE, hier.names = FALSE, use.dnns = FALSE) {
    if (is.character(hier.names))
        hier.names <- match.arg(hier.names, c("before", "after",
                                              "none"))
    else hier.names <- if (hier.names)
        "before"
    else "no"
    arg.list <- list(...)
    if (is.list(arg.list[[1]]) && !is.data.frame(arg.list[[1]])) {
        if (length(arg.list) != 1)
            stop("can only supply one list-valued argument for ...")
        if (make.names)
            stop("cannot have make.names=TRUE with a list argument")
        arg.list <- arg.list[[1]]
        have.list.arg <- TRUE
    }
    else {
        N <- max(1, sapply(list(...), function(x) length(dim(x))))
        have.list.arg <- FALSE
    }
    if (any(discard <- sapply(arg.list, is.null)))
        arg.list <- arg.list[!discard]
    if (length(arg.list) == 0)
        return(NULL)
    N <- max(1, sapply(arg.list, function(x) length(dim(x))))
    if (!is.null(rev.along))
        along <- N + 1 - rev.along
    if (along < 1 || along > N || (along > floor(along) && along <
                                   ceiling(along))) {
        N <- N + 1
        along <- max(1, min(N + 1, ceiling(along)))
    }
    if (length(along) > 1 || along < 1 || along > N + 1)
        stop(paste("\"along\" must specify one dimension of the array,",
                   "or interpolate between two dimensions of the array",
                   sep = "\n"))
    if (!force.array && N == 2) {
        if (!have.list.arg) {
            if (along == 2)
                return(cbind(...))
            if (along == 1)
                return(rbind(...))
        }
        else {
            if (along == 2)
                return(do.call("cbind", arg.list))
            if (along == 1)
                return(do.call("rbind", arg.list))
        }
    }
    if (along > N || along < 0)
        stop("along must be between 0 and ", N)
    pre <- seq(from = 1, len = along - 1)
    post <- seq(to = N - 1, len = N - along)
    perm <- c(seq(len = N)[-along], along)
    arg.names <- names(arg.list)
    if (is.null(arg.names))
        arg.names <- rep("", length(arg.list))
    if (is.character(new.names)) {
        arg.names[seq(along = new.names)[nchar(new.names) > 0]] <- new.names[nchar(new.names) >
                                                                                 0]
        new.names <- NULL
    }
    if (any(arg.names == "")) {
        if (make.names) {
            dot.args <- match.call(expand.dots = FALSE)$...
            if (is.call(dot.args) && identical(dot.args[[1]],
                                               as.name("list")))
                dot.args <- dot.args[-1]
            arg.alt.names <- arg.names
            for (i in seq(along = arg.names)) {
                if (arg.alt.names[i] == "") {
                    if (object.size(dot.args[[i]]) < 1000) {
                        arg.alt.names[i] <- paste(deparse(dot.args[[i]],
                                                          40), collapse = ";")
                    }
                    else {
                        arg.alt.names[i] <- paste("X", i, sep = "")
                    }
                    arg.names[i] <- arg.alt.names[i]
                }
            }
        }
        else {
            arg.alt.names <- arg.names
            arg.alt.names[arg.names == ""] <- paste("X", seq(along = arg.names),
                                                    sep = "")[arg.names == ""]
        }
    }
    else {
        arg.alt.names <- arg.names
    }
    use.along.names <- any(arg.names != "")
    names(arg.list) <- arg.names
    arg.dimnames <- matrix(vector("list", N * length(arg.names)),
                           nrow = N, ncol = length(arg.names))
    dimnames(arg.dimnames) <- list(NULL, arg.names)
    arg.dnns <- matrix(vector("list", N * length(arg.names)),
                       nrow = N, ncol = length(arg.names))
    dimnames(arg.dnns) <- list(NULL, arg.names)
    dimnames.new <- vector("list", N)
    arg.dim <- matrix(integer(1), nrow = N, ncol = length(arg.names))
    for (i in seq(len = length(arg.list))) {
        m <- arg.list[[i]]
        m.changed <- FALSE
        if (is.data.frame(m)) {
            m <- as.matrix(m)
            m.changed <- TRUE
        }
        else if (!is.array(m) && !is.null(m)) {
            if (!is.atomic(m))
                stop("arg '", arg.alt.names[i], "' is non-atomic")
            dn <- names(m)
            m <- as.array(m)
            if (length(dim(m)) == 1 && !is.null(dn))
                dimnames(m) <- list(dn)
            m.changed <- TRUE
        }
        new.dim <- dim(m)
        if (length(new.dim) == N) {
            if (!is.null(dimnames(m))) {
                arg.dimnames[, i] <- dimnames(m)
                if (use.dnns && !is.null(names(dimnames(m))))
                    arg.dnns[, i] <- as.list(names(dimnames(m)))
            }
            arg.dim[, i] <- new.dim
        }
        else if (length(new.dim) == N - 1) {
            if (!is.null(dimnames(m))) {
                arg.dimnames[-along, i] <- dimnames(m)
                if (use.dnns && !is.null(names(dimnames(m))))
                    arg.dnns[-along, i] <- as.list(names(dimnames(m)))
                dimnames(m) <- NULL
            }
            arg.dim[, i] <- c(new.dim[pre], 1, new.dim[post])
            if (any(perm != seq(along = perm))) {
                dim(m) <- c(new.dim[pre], 1, new.dim[post])
                m.changed <- TRUE
            }
        }
        else {
            stop("'", arg.alt.names[i], "' does not fit: should have `length(dim())'=",
                 N, " or ", N - 1)
        }
        if (any(perm != seq(along = perm)))
            arg.list[[i]] <- aperm(m, perm)
        else if (m.changed)
            arg.list[[i]] <- m
    }
    conform.dim <- arg.dim[, 1]
    for (i in seq(len = ncol(arg.dim))) {
        if (any((conform.dim != arg.dim[, i])[-along])) {
            stop("arg '", arg.alt.names[i], "' has dims=", paste(arg.dim[,
                                                                         i], collapse = ", "), "; but need dims=", paste(replace(conform.dim,
                                                                                                                                 along, "X"), collapse = ", "))
        }
    }
    if (N > 1)
        for (dd in seq(len = N)[-along]) {
            for (i in (if (use.first.dimnames)
                seq(along = arg.names)
                else rev(seq(along = arg.names)))) {
                if (length(arg.dimnames[[dd, i]]) > 0) {
                    dimnames.new[[dd]] <- arg.dimnames[[dd, i]]
                    if (use.dnns && !is.null(arg.dnns[[dd, i]]))
                        names(dimnames.new)[dd] <- arg.dnns[[dd,
                                                             i]]
                    break
                }
            }
        }
    for (i in seq(len = length(arg.names))) {
        if (arg.dim[along, i] > 0) {
            dnm.along <- arg.dimnames[[along, i]]
            if (length(dnm.along) == arg.dim[along, i]) {
                use.along.names <- TRUE
                if (hier.names == "before" && arg.names[i] !=
                    "")
                    dnm.along <- paste(arg.names[i], dnm.along,
                                       sep = ".")
                else if (hier.names == "after" && arg.names[i] !=
                         "")
                    dnm.along <- paste(dnm.along, arg.names[i],
                                       sep = ".")
            }
            else {
                if (arg.dim[along, i] == 1)
                    dnm.along <- arg.names[i]
                else if (arg.names[i] == "")
                    dnm.along <- rep("", arg.dim[along, i])
                else dnm.along <- paste(arg.names[i], seq(length = arg.dim[along,
                                                                           i]), sep = "")
            }
            dimnames.new[[along]] <- c(dimnames.new[[along]],
                                       dnm.along)
        }
        if (use.dnns) {
            dnn <- unlist(arg.dnns[along, ])
            if (length(dnn)) {
                if (!use.first.dimnames)
                    dnn <- rev(dnn)
                names(dimnames.new)[along] <- dnn[1]
            }
        }
    }
    if (!use.along.names)
        dimnames.new[along] <- list(NULL)
    out <- array(unlist(arg.list, use.names = FALSE),
                 dim = c(arg.dim[-along, 1], sum(arg.dim[along, ])),
                 dimnames = dimnames.new[perm])
    if (any(order(perm) != seq(along = perm)))
        out <- aperm(out, order(perm))
    if (!is.null(new.names) && is.list(new.names)) {
        for (dd in seq(len = N)) {
            if (!is.null(new.names[[dd]])) {
                if (length(new.names[[dd]]) == dim(out)[dd])
                    dimnames(out)[[dd]] <- new.names[[dd]]
                else if (length(new.names[[dd]]))
                    warning(paste("Component ", dd, " of new.names ignored: has length ",
                                  length(new.names[[dd]]), ", should be ",
                                  dim(out)[dd], sep = ""))
            }
            if (use.dnns && !is.null(names(new.names)) && names(new.names)[dd] !=
                "")
                names(dimnames(out))[dd] <- names(new.names)[dd]
        }
    }
    if (use.dnns && !is.null(names(dimnames(out))) && any(i <- is.na(names(dimnames(out)))))
        names(dimnames(out))[i] <- ""
    out
}

