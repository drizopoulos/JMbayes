anova.JMbayes <-
function (object, ...) {
    ancall <- sys.call()
    names <- unlist(lapply(as.list(ancall[-1L]), deparse))
    if (length(objects <- list(object, ...)) == 1L)
        stop("anova is only used to compare multiple joint models.\n")
    ns <- sapply(objects, function (obj) length(obj$y$Time))
    if (any(ns != ns[1L])) 
        stop("models were not all fitted to the same size of dataset.")
    table <- as.data.frame(t(sapply(objects, function (obj) c(df = length(unlist(obj$postMeans)), LPML = obj$LPML, 
                                                DIC = obj$DIC, pD = obj$pD))))
    row.names(table) <- names
    table
}
