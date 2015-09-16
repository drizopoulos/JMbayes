makepredictcall.dbs <-
function (var, call) {
    if (as.character(call)[1L] != "dbs") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept", "eps")]
    xxx <- call[1L:2L]
    xxx[names(at)] <- at
    xxx
}
