format.perc <-
function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, 
                 scientific = FALSE, digits = digits), "%")
}
