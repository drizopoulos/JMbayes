plot.survfitJM2 <- function (x, estimator = c("both", "mean", "median"), 
          which = NULL, fun = NULL, invlink = NULL, conf.int = FALSE, 
          fill.area = FALSE, col.area = "grey", col.abline = "black", 
          col.points = "black", add.last.time.axis.tick = FALSE, 
          include.y = FALSE, main = NULL, xlab = NULL, ylab = NULL, 
          ylab2 = NULL, lty = NULL, col = NULL, lwd = NULL, pch = NULL, 
          ask = NULL, legend = FALSE, ..., cex.axis.z = 1, cex.lab.z = 1, 
          xlim = NULL, ylim2 = NULL) 
{
  estimator <- match.arg(estimator)
  fun <- if (!is.null(fun)) 
    match.fun(fun)
  if (is.null(which)) 
    which <- seq_along(x$summaries)
  if (conf.int && is.null(x$success.rate)) {
    warning("a confidence interval can be included only when argument", 
            "'simulate' of survfitJMBUGS() was set to TRUE.")
    conf.int <- FALSE
  }
  if (is.null(ask)) 
    ask <- prod(par("mfcol")) < length(which)
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (is.null(main)) {
    main <- paste("Subject", names(x$summaries))
    names(main) <- names(x$summaries)
  }
  if (is.null(xlab)) 
    xlab <- rep("Time", length(which))
  if (is.null(ylab)) {
    ylab <- if (is.null(fun)) 
      if (!include.y) 
        rep(expression(paste("Pr(", T[i] >= u, 
                             " | ", T[i] > t, ", ", tilde(y)[i](t), 
                             ")", sep = " ")), length(which))
    else rep("Event-free Probability", length(which))
    else rep("", length(which))
  }
  if (is.null(ylab2)) 
    ylab2 <- x$nameY
  if (!is.null(x$success.rate)) {
    if (is.null(col)) 
      col <- switch(estimator, both = c(2, 3, 1, 1), mean = c(2, 
                                                              1, 1), median = c(3, 1, 1))
    if (is.null(lty)) 
      lty <- switch(estimator, both = c(1, 1, 2, 2), mean = c(1, 
                                                              2, 2), median = c(1, 2, 2))
    if (is.null(lwd)) 
      lwd <- switch(estimator, both = c(1, 1, 1, 1), mean = c(1, 
                                                              1, 1), median = c(1, 1, 1))
  }
  else {
    col <- lty <- lwd <- 1
  }
  if (is.null(pch)) 
    pch <- 8
  for (i in seq_along(which)) {
    ii <- which[i]
    r <- x$summaries[[ii]]
    r <- if (!is.null(x$success.rate)) {
      rbind(cbind(c(0, x$last.time[ii]), matrix(1, 2, 4)), 
            r)
    }
    else {
      rbind(cbind(c(0, x$last.time[ii]), matrix(1, 2, 1)), 
            r)
    }
    if (!is.null(fun) && is.function(fun)) 
      r[, 2:ncol(r)] <- fun(r[, 2:ncol(r)])
    if (!is.null(x$success.rate) && estimator == "mean") 
      r <- r[, -3]
    if (!is.null(x$success.rate) && estimator == "median") 
      r <- r[, -2]
    if (!conf.int && !is.null(x$success.rate)) {
      exc <- c(ncol(r) - 1, ncol(r))
      r <- r[, -exc, drop = FALSE]
      col <- col[-exc]
      lty <- lty[-exc]
      lwd <- lwd[-exc]
    }
    ylim <- if (is.null(fun)) 
      c(0, 1)
    else {
      rr <- r[, -1, drop = FALSE]
      range(rr[is.finite(rr)])
    }
    if (!include.y) {
      matplot(r[, 1], r[, -1, drop = FALSE], type = "l", 
              col = col, lwd = lwd, lty = lty, ylim = ylim, 
              main = main[ii], xlab = xlab[i], ylab = ylab[i], 
              ...)
      if (fill.area) {
        polygon(c(r[, 1], rev(r[, 1])), c(r[, ncol(r) - 
                                              1], rev(r[, ncol(r)])), col = col.area, border = "transparent")
        matlines(r[, 1], r[, -1, drop = FALSE], type = "l", 
                 col = col, lwd = lwd, lty = lty)
      }
    }
    else {
      if(!is.null(invlink)) 
        x$y[[ii]] <- invlink(x$y[[ii]]) 
      else x$y[[ii]] <- x$y[[ii]]
      oldmar <- par("mar")
      par(mar = c(5, 4, 5, 4))
      lt <- x$last.time[ii]
      r. <- r[r[, 1] >= lt, ]
      rng <- if (is.null(xlim)) 
        range(x$obs.times[[ii]], x$survTimes)
      else xlim
      rng2 <- if(is.null(ylim2))
        x$ry
      else ylim2 
        
      plot(x$obs.times[[ii]], x$y[[ii]], xlim = rng, ylim = rng2, 
           xlab = xlab[i], ylab = ylab2, pch = pch, col = col.points, 
           ...)
      
      ff <- if (!is.null(invlink)) 
        invlink(x$fitted.y[[ii]]) #& invlink(x$y[[ii]])
      else x$fitted.y[[ii]] #& x$y[[ii]]
      lines(x$fitted.times[[ii]], ff, col = col, lwd = lwd)
      abline(v = lt, lty = 3, col = col.abline)
      par(new = TRUE)
      matplot(r.[, 1], r.[, -1, drop = FALSE], type = "l", 
              col = col, lwd = lwd, lty = lty, ylim = ylim, 
              main = main[ii], xlim = rng, ylab = "", 
              xlab = "", axes = FALSE, yaxs = "i", 
              ...)
      axis(4, las = 2, cex.axis = cex.axis.z)
      if (fill.area) {
        polygon(c(r.[, 1], rev(r.[, 1])), c(r.[, ncol(r.) - 
                                                 1], rev(r.[, ncol(r.)])), col = col.area, border = "transparent")
        matlines(r.[, 1], r.[, -1, drop = FALSE], type = "l", 
                 col = col, lwd = lwd, lty = lty)
      }
      mtext(ylab[i], 4, 2, cex = cex.lab.z)
      par(mar = oldmar, new = FALSE)
    }
    if (add.last.time.axis.tick) 
      axis(1, at = round(x$last.time[ii], 1))
    if (legend) {
      lab <- switch(estimator, both = c("Mean", "Median"), 
                    mean = "Mean", median = "Median")
      if (conf.int) 
        lab <- c(lab, "Lower limit", "Upper limit")
      legend("left", lab, lwd = lwd, lty = lty, col = col, 
             bty = "n", ...)
    }
  }
  invisible()
}
