runDynPred <- function (type = c("JM", "lme")) {
    type <- match.arg(type)
    if (requireNamespace("shiny", quietly = TRUE)) {
        if (type == "JM") {
            shiny::runApp(system.file("shiny_app_JM", package = "JMbayes"))
        } else {
            shiny::runApp(system.file("shiny_app_lme", package = "JMbayes"))
        }
    } else {
        cat("\npackage 'shiny' is not installed...")
    }
}
