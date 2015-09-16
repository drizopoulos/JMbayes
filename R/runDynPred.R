runDynPred <-
function () {
    if (requireNamespace("shiny", quietly = TRUE)) {
        shiny::runApp(system.file("Demo", package = "JMbayes"))
    } else {
        cat("\npackage 'shiny' is not installed...")
    }
}
