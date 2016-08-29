runDynPred <-
function () {
    if (requireNamespace("shiny", quietly = TRUE)) {
        shiny::runApp(system.file("demo", package = "JMbayes"))
    } else {
        cat("\npackage 'shiny' is not installed...")
    }
}
