makepredictcall.ibs <-
function (var, call) {
    if (as.character(call)[1L] != "ibs") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept", "from", "weight.fun")]
    xxx <- call[1L:2L]
    xxx[names(at)] <- at
    xxx
}
