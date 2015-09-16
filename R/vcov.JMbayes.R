vcov.JMbayes <-
function (object, ...) {
    vmat <- object$vcov
    (vmat + t(vmat)) / 2
}
