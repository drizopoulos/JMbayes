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
        timeTerms <- if (length(terms.labs_Z)) grep(terms.labs_Z[1L], colnames(X), fixed = TRUE)
        which_td <- unname(which(apply(X, 2, check_td, id = id)))
        all_TDterms <- unique(c(timeTerms, which_td))
        baseline <- seq_len(ncol(X))[-all_TDterms]
        ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions,
                                               nams2 = colnames(X)))
        ind_colmns2 <- seq_len(ncol(X))
        ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
        data.id <- data[!duplicated(id), ]
        Xhc <- if (length(terms.labs_Z)) {
            data.id[[terms.labs_Z[1L]]] <- 1
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

dwish <- function (W, S, v, log = FALSE) {
    if (!is.matrix(S))
        S <- matrix(S)
    if (!is.matrix(W))
        W <- matrix(W)
    k <- nrow(S)
    log.gammapart <- 0
    for (i in seq_len(k)) {
        log.gammapart <- log.gammapart + lgamma(0.5 * (v + 1 - i))
    }
    log.denom <- log.gammapart + (0.5 * v * k) * log(2) + (0.25 * k * (k - 1)) * log(pi)
    log.detS <- c(determinant(S, logarithm = TRUE)$modulus)
    log.detW <- c(determinant(W, logarithm = TRUE)$modulus)
    trace <- sum(diag(solve(S, W)))
    log.num <- - 0.5 * (v * log.detS - (v - k - 1) * log.detW + trace)
    if (log)
        log.num - log.denom
    else
        exp(log.num - log.denom)
}

gaussKronrod <- function (k = 15L) {
    sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397,
            0, 0.405845151377397, 0.741531185599394, 0.949107912342758,
            -0.991455371120813, -0.864864423359769, -0.586087235467691,
            -0.207784955007898, 0.207784955007898, 0.586087235467691,
            0.864864423359769, 0.991455371120813)
    wk15 <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785,
              0.209482141084728, 0.190350578064785, 0.140653259715526,
              0.0630920926299786, 0.0229353220105292, 0.10479001032225,
              0.169004726639268, 0.204432940075299, 0.204432940075299,
              0.169004726639268, 0.10479001032225, 0.0229353220105292)
    wk7 <- c(0.12948496616887, 0.279705391489277, 0.381830050505119,
             0.417959183673469, 0.381830050505119, 0.279705391489277,
             0.12948496616887)
    out <- if (k == 7L) {
        list(sk = sk[1:7], wk = wk7)
    }
    else {
        list(sk = sk, wk = wk15)
    }
    out
}

gaussLegendre <- function (k = 15) {
    if (!k %in% c(15, 16, 17, 18, 19, 20, 32))
        stop("k should be 15:20 or 32.")
    out <- if (k == 15) {
        list(sk = c(0, 0.201194093997435, -0.201194093997435,
                    0.394151347077563, -0.394151347077563, 0.570972172608539,
                    -0.570972172608539, 0.72441773136017, -0.72441773136017,
                    0.848206583410427, -0.848206583410427, 0.937273392400706,
                    -0.937273392400706, 0.987992518020485, -0.987992518020485),
             wk = c(0.202578241925561, 0.198431485327112, 0.198431485327112,
                    0.186161000015562, 0.186161000015562, 0.166269205816994,
                    0.166269205816994, 0.139570677926154, 0.139570677926154,
                    0.107159220467172, 0.107159220467172, 0.0703660474881081,
                    0.0703660474881081, 0.0307532419961173, 0.0307532419961173))
    }
    else if (k == 16) {
        list(sk = c(0.0950125098376374, -0.0950125098376374,
                    0.281603550779259, -0.281603550779259, 0.458016777657227,
                    -0.458016777657227, 0.617876244402644, -0.617876244402644,
                    0.755404408355003, -0.755404408355003, 0.865631202387832,
                    -0.865631202387832, 0.944575023073233, -0.944575023073233,
                    0.98940093499165, -0.98940093499165), wk = c(0.189450610455069,
                                                                 0.189450610455069, 0.182603415044924, 0.182603415044924,
                                                                 0.169156519395003, 0.169156519395003, 0.149595988816577,
                                                                 0.149595988816577, 0.124628971255534, 0.124628971255534,
                                                                 0.0951585116824928, 0.0951585116824928, 0.0622535239386479,
                                                                 0.0622535239386479, 0.0271524594117541, 0.0271524594117541))
    }
    else if (k == 17) {
        list(sk = c(0, 0.178484181495848, -0.178484181495848,
                    0.351231763453876, -0.351231763453876, 0.512690537086477,
                    -0.512690537086477, 0.657671159216691, -0.657671159216691,
                    0.781514003896801, -0.781514003896801, 0.880239153726986,
                    -0.880239153726986, 0.950675521768768, -0.950675521768768,
                    0.990575475314417, -0.990575475314417), wk = c(0.179446470356207,
                                                                   0.176562705366993, 0.176562705366993, 0.16800410215645,
                                                                   0.16800410215645, 0.15404576107681, 0.15404576107681,
                                                                   0.135136368468525, 0.135136368468525, 0.111883847193404,
                                                                   0.111883847193404, 0.0850361483171792, 0.0850361483171792,
                                                                   0.0554595293739872, 0.0554595293739872, 0.0241483028685479,
                                                                   0.0241483028685479))
    }
    else if (k == 18) {
        list(sk = c(0.0847750130417353, -0.0847750130417353,
                    0.251886225691505, -0.251886225691505, 0.411751161462843,
                    -0.411751161462843, 0.559770831073948, -0.559770831073948,
                    0.691687043060353, -0.691687043060353, 0.803704958972523,
                    -0.803704958972523, 0.892602466497556, -0.892602466497556,
                    0.955823949571398, -0.955823949571398, 0.991565168420931,
                    -0.991565168420931), wk = c(0.169142382963144, 0.169142382963144,
                                                0.164276483745833, 0.164276483745833, 0.154684675126265,
                                                0.154684675126265, 0.140642914670651, 0.140642914670651,
                                                0.122555206711478, 0.122555206711478, 0.100942044106287,
                                                0.100942044106287, 0.0764257302548891, 0.0764257302548891,
                                                0.0497145488949698, 0.0497145488949698, 0.0216160135264833,
                                                0.0216160135264833))
    }
    else if (k == 19) {
        list(sk = c(0, 0.160358645640225, -0.160358645640225,
                    0.31656409996363, -0.31656409996363, 0.464570741375961,
                    -0.464570741375961, 0.600545304661681, -0.600545304661681,
                    0.720966177335229, -0.720966177335229, 0.822714656537143,
                    -0.822714656537143, 0.903155903614818, -0.903155903614818,
                    0.96020815213483, -0.96020815213483, 0.992406843843584,
                    -0.992406843843584), wk = c(0.161054449848784, 0.158968843393954,
                                                0.158968843393954, 0.15276604206586, 0.15276604206586,
                                                0.142606702173607, 0.142606702173607, 0.128753962539336,
                                                0.128753962539336, 0.111566645547334, 0.111566645547334,
                                                0.09149002162245, 0.09149002162245, 0.0690445427376412,
                                                0.0690445427376412, 0.0448142267656996, 0.0448142267656996,
                                                0.0194617882297265, 0.0194617882297265))
    }
    else if (k == 20) {
        list(sk = c(0.0765265211334973, -0.0765265211334973,
                    0.227785851141645, -0.227785851141645, 0.37370608871542,
                    -0.37370608871542, 0.510867001950827, -0.510867001950827,
                    0.636053680726515, -0.636053680726515, 0.746331906460151,
                    -0.746331906460151, 0.839116971822219, -0.839116971822219,
                    0.912234428251326, -0.912234428251326, 0.963971927277914,
                    -0.963971927277914, 0.993128599185095, -0.993128599185095),
             wk = c(0.152753387130726, 0.152753387130726, 0.149172986472604,
                    0.149172986472604, 0.142096109318382, 0.142096109318382,
                    0.131688638449177, 0.131688638449177, 0.118194531961518,
                    0.118194531961518, 0.10193011981724, 0.10193011981724,
                    0.0832767415767048, 0.0832767415767048, 0.0626720483341091,
                    0.0626720483341091, 0.0406014298003869, 0.0406014298003869,
                    0.0176140071391521, 0.0176140071391521))
    }
    else {
        list(sk = c(0.0483076656877383, -0.0483076656877383,
                    0.144471961582796, -0.144471961582796, 0.239287362252137,
                    -0.239287362252137, 0.331868602282128, -0.331868602282128,
                    0.421351276130635, -0.421351276130635, 0.506899908932229,
                    -0.506899908932229, 0.587715757240762, -0.587715757240762,
                    0.663044266930215, -0.663044266930215, 0.73218211874029,
                    -0.73218211874029, 0.794483795967942, -0.794483795967942,
                    0.84936761373257, -0.84936761373257, 0.896321155766052,
                    -0.896321155766052, 0.93490607593774, -0.93490607593774,
                    0.964762255587506, -0.964762255587506, 0.985611511545268,
                    -0.985611511545268, 0.997263861849482, -0.997263861849482),
             wk = c(0.0965400885147278, 0.0965400885147278, 0.0956387200792749,
                    0.0956387200792749, 0.0938443990808046, 0.0938443990808046,
                    0.0911738786957639, 0.0911738786957639, 0.0876520930044038,
                    0.0876520930044038, 0.0833119242269467, 0.0833119242269467,
                    0.0781938957870703, 0.0781938957870703, 0.0723457941088485,
                    0.0723457941088485, 0.0658222227763618, 0.0658222227763618,
                    0.0586840934785355, 0.0586840934785355, 0.0509980592623762,
                    0.0509980592623762, 0.0428358980222267, 0.0428358980222267,
                    0.0342738629130214, 0.0342738629130214, 0.0253920653092621,
                    0.0253920653092621, 0.0162743947309057, 0.0162743947309057,
                    0.0070186100094701, 0.0070186100094701))
    }
    ord <- order(out$sk)
    f <- function(x, ord) x[ord]
    lapply(out, f, ord = ord)
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

