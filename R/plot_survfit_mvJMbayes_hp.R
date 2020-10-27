##################################################
## Project: PhD - to add to 
## https://github.com/drizopoulos/JMbayes/blob/3ba136bd4f04cc53554196ca3a48dd06793e95bb/R/survfitJM.mvJMbayes.R
## Script purpose: amending github and hardcode over it 
## Date: 2020-03-31
## Author: Harry Parr
##################################################


plot_survfit_mvJMbayes_hp <- function (x, split = c(1, 1), splitr = c(1,1), which_subjects = NULL, which_outcomes = NULL,
                                    surv_in_all = TRUE, include.y = TRUE, fun = NULL,
                                    abline = NULL, theinvlink = T, yexp = FALSE, 
                                    main = NULL, xlab = "Time", ylab = NULL, 
                                    zlab = "Event-Free Probability",
                                    include_CI = TRUE, fill_area_CI = TRUE, 
                                    col_points = "black", pch_points = 1,
                                    col_lines = "red", col_lines_CI = "black", 
                                    col_fill_CI = "lightgrey",
                                    lwd_lines = 2, lty_lines_CI = 2, 
                                    cex_xlab = 1, cex_ylab = 1, cex_zlab = 1, cex_main = 1,
                                    cex_axis = 1, ...) {
  #families <- x$families
  families <- x$families
  #families <- if (theinvlink) x$families <- exp else (x$families)
  #families <- if (!is.null(theinvlink)) x$families <- theinvlink else (x$families)
  ylim <- NULL
  respVars <- x$respVars
  N <- length(x$y)
  n_outcomes <- length(x$y[[1]])
  xlim <- range(unlist(x$obs.times), x$survTimes)
  add_surv <- function (xaxis = TRUE, outer = TRUE) {
    summ <- x$summaries[[i]]
    times <- c(x$last.time[[i]], summ[, "times"])
    surv <- c(1, summ[, "Mean"])
    low <- c(1, summ[, "Lower"])
    upp <- c(1, summ[, "Upper"])
    if (!is.null(fun) && is.function(fun)) {
      surv <- fun(surv)
      low <- fun(low)
      upp <- fun(upp)
    }
    plot(times, surv, col = col_lines, type = "l", ylim = c(0.0, 1),
         lwd = 2, xlim = xlim, axes = FALSE)
    box()
    if (include_CI && fill_area_CI) {
      polygon(c(times, rev(times)), c(low, rev(upp)), border = NA, col = col_fill_CI)
    }
    if (include_CI) {
      lines(times, low, lwd = lwd_lines, lty = lty_lines_CI, col = col_lines_CI)
      lines(times, upp, lwd = lwd_lines, lty = lty_lines_CI, col = col_lines_CI)
    }
    lines(times, surv, lwd = lwd_lines, col = col_lines)
    abline(v = x$last.time[[i]], lty = 2)
    if (!is.null(abline)) {
      abline(v = abline$v, lty = abline$lty, lwd = abline$lwd, col = abline$col)
      if (xaxis) axis(1, at = round(abline$v, 1), cex = cex_axis)
    }
    if (xaxis) axis(1, cex.axis = cex_axis)
    axis(4, cex.axis = cex_axis)
    mtext(zlab, side = 4, line = 1.8, outer = outer, cex = cex_zlab)
  }
  if (is.null(ylab))
    ylab <- respVars
  if (is.null(main))
    main <- paste("Subject", names(x$y))
  valid_subjects <- seq_len(N)
  if (!is.null(which_subjects)) {
    if (!all(which_subjects %in% valid_subjects)) {
      stop("'which_subjects' must be an integer vector with possible values: ",
           paste(valid_subjects, collapse = ", "))
    } else {
      valid_subjects <- which_subjects
    }
  }
  valid_outcomes <- seq_len(n_outcomes)
  if (!is.null(which_outcomes)) {
    if (!all(which_outcomes %in% valid_outcomes)) {
      stop("'which_outcomes' must be an integer vector with possible values: ",
           paste(valid_outcomes, collapse = ", "))
    } else {
      valid_outcomes <- which_outcomes
    }
  }
  for (i in valid_subjects) {
    opar <- par(no.readonly = TRUE, mfcol = split, mfrow = splitr, oma = c(3, 3, 2, 3), 
                mar = c(0, 0, 0, 0), mgp = c(3, 0.4, 0), tcl = -0.25)
    if (include.y) {
      for (j in valid_outcomes) {
        obs_times <- x$obs.times[[i]]
         # y <- x$y[[i]][[j]]
        y <- if (yexp) sapply(x$y[[i]][[j]] , exp) else x$y[[i]][[j]]
        #y <- if(!is.null(theinvlink)) theinvlink(x$y[[i]][[j]]) else x$y[[i]][[j]]
        if (fact_y <- is.factor(y)) {
          lvy <- levels(y)
          y <- as.numeric(y == levels(y)[2L])
          ylim <- c(-0.2, 1.2)
        }
        fitted_times <- x$fitted.times[[i]]
         # fitted_y <- families[[j]]$linkinv(x$fitted.y[[i]][[j]])
        fitted_y <- if(yexp) sapply(families[[j]]$linkinv(x$fitted.y[[i]][[j]]), exp) else families[[j]]$linkinv(x$fitted.y[[i]][[j]])
         
         #families <- if (theinvlink) x$families <- exp else (x$families)
         #families <- if (!is.null(theinvlink)) x$families <- theinvlink else (x$families)
        #fitted_y <- if(!is.null(theinvlink)) families[[j]]$theinvlink(x$fitted.y[[i]][[j]]) else families[[j]]$linkinv(x$fitted.y[[i]][[j]])
        #if(splitr != c(1,1)) par(mfrow=splitr) else par(mfrow=c(1,1))
        par(opar)
        plot(obs_times, y, ylim = ylim, ylab = respVars[i], xlim = xlim, axes = FALSE,
             type = "n")
        box()
        points(obs_times, y, col = col_points, pch = pch_points)
        if (fact_y) axis(2, at = 0:1, labels = lvy, cex.axis = cex_axis) else axis(2, cex.axis = cex_axis)
        if (add_xaxis <- par()$mfcol[1L] == j) axis(1, cex.axis = cex_axis)
        lines(fitted_times, fitted_y, lwd = lwd_lines, col = col_lines)
        abline(v = x$last.time[[i]], lty = 2)
        mtext(ylab[j], side = 2, line = 1.8, cex = cex_ylab)
        ylim <- NULL
        if (surv_in_all) {
          par(new = TRUE)
          add_surv(add_xaxis)
        }
        if (added_xlab <- all(par()$mfcol == c(1, 1))) {
          mtext(xlab, side = 1, line = 1.5, outer = TRUE, cex = cex_xlab)
          mtext(main[i], side = 3, line = 0.8, outer = TRUE, cex = cex_main)
        }
      }
    }
    if (!include.y || !surv_in_all) {
      add_surv(outer = FALSE)
    }
    if ((!include.y || !added_xlab)) {
      mtext(xlab, side = 1,line = 1.5, outer = TRUE, cex = cex_xlab)
      mtext(main[i], side = 3, line = 0.8, outer = TRUE, cex = cex_main)
    }
    if (length(valid_outcomes) == 1) {
      axis(1, cex.axis = cex_axis)
      if (!is.null(abline)) 
        axis(1, at = round(abline$v, 1), cex = cex_axis)
    }
    par(opar)
  }
  invisible()
}

##################################################
##################################################
## Section: Testing plot
##################################################
#mvJMFit_tveffectN3071$model_info$mvglmer_components$data
# ND <- mvJMFit_tveffectN3071$model_info$mvglmer_components$data[mvJMFit_tveffectN3071$model_info$mvglmer_components$data$trialno %in% c(2, 25, 81), ]
# ND_182024_172024
# #sprobs <- survfitJM(mvJMFit_tveffectN3071, ND_182024_172024)
# #class(sprobs) <- "survfit.JMbayes"
# plot(sprobs)
# plot_survfit_mvJMbayes_hp(sprobs)
# plot_survfit_mvJMbayes_hp(sprobs, theinvlink = T)
# plot_survfit_mvJMbayes_hp(sprobs, theinvlink = T, yexp=TRUE, ylab="PSA", main="")


plot_survfit_mvJMbayes_hp(sprobs, split = c(2, 1), surv_in_all = TRUE, xlab="Time (years)", main="",
                          ylab="PSA ng/mL", yexp = TRUE, which_subjects = c(2,1), splitr = c(2,1))


# N <- nrow(ND_182024_172024)
# dyn_sprobs <- vector("list", N)
# dyn_sprobsA <- vector("list", 16)
# dyn_sprobsB <- vector("list", 16)
# 
# par(mfrow=c(2,1))
# ND_172024 <- subset(ND_182024_172024, trialno == 172024)
# for (i in seq_len(16)) {
#   par(mfrow=c(2,1))
#   dyn_sprobsA[[i]] <- survfitJM(mvJMFit_tveffectN3071, ND_182024[1:i, ], 
#                                survTimes = seq(0, 10, length.out = 85))
#   #plot(dyn_sprobs[[i]], split = c(1, 1), surv_in_all = TRUE)
#   plot_survfit_mvJMbayes_hp(dyn_sprobsA[[i]], split = c(2, 1), surv_in_all = TRUE, xlab="Time (years)", main="Patient A",
#                             ylab="PSA ng/mL", yexp = TRUE, which_subjects = 1)
#   dyn_sprobsB[[i]] <- survfitJM(mvJMFit_tveffectN3071, ND_172024, 
#                                survTimes = seq(0, 10, length.out = 85))
#   plot_survfit_mvJMbayes_hp(dyn_sprobsB[[i]], split = c(2, 1), surv_in_all = TRUE, xlab="Time (years)", main="Patient B",
#                             ylab="PSA ng/mL", yexp = TRUE, which_subjects = 1)
# }

##################################################
## Section: GIFs for DPs
##################################################

# animation::saveGIF(
# for (i in seq_len(16)) {
#   par(mfrow=c(2,1))
#   dyn_sprobsA[[i]] <- survfitJM(mvJMFit_tveffectN3071, ND_182024[1:i, ], 
#                                 survTimes = seq(0, 10, length.out = 85))
#   #plot(dyn_sprobs[[i]], split = c(1, 1), surv_in_all = TRUE)
#   plot_survfit_mvJMbayes_hp(dyn_sprobsA[[i]], split = c(2, 1), surv_in_all = TRUE, xlab="Time (years)", main="Patient A",
#                             ylab="PSA ng/mL", yexp = TRUE, which_subjects = 1)
#   dyn_sprobsB[[i]] <- survfitJM(mvJMFit_tveffectN3071, ND_172024, 
#                                 survTimes = seq(0, 10, length.out = 85))
#   plot_survfit_mvJMbayes_hp(dyn_sprobsB[[i]], split = c(2, 1), surv_in_all = TRUE, xlab="Time (years)", main="Patient B",
#                             ylab="PSA ng/mL", yexp = TRUE, which_subjects = 1)
# }
# , "dynamic_preds_mvJM.gif", ani.width=1000, ani.height=600)
