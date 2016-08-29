\name{plot.survfit.JMbayes}

\alias{plot.survfit.JMbayes}

\title{Plot Method for survfit.JMbayes Objects}

\description{
  Produces plots of conditional probabilities of survival.
}

\usage{
\method{plot}{survfit.JMbayes}(x, estimator = c("both", "mean", "median"), 
    which = NULL, fun = NULL, invlink = NULL, conf.int = FALSE, 
    fill.area = FALSE, col.area = "grey", col.abline = "black", col.points = "black",
    add.last.time.axis.tick = FALSE, include.y = FALSE, main = NULL, 
    xlab = NULL, ylab = NULL, ylab2 = NULL, lty = NULL, col = NULL, 
    lwd = NULL, pch = NULL, ask = NULL, legend = FALSE, \dots,
    cex.axis.z = 1, cex.lab.z = 1, xlim = NULL)
}

\arguments{
  \item{x}{an object inheriting from class \code{survfit.JMbayes}.}
  \item{estimator}{character string specifying, whether to include in the plot the mean of the conditional probabilities of survival,
    the median or both. The mean and median are taken as estimates of these conditional probabilities over the M replications of the
    Monte Carlo scheme described in \code{\link{survfitJM}}.}
  \item{which}{a numeric or character vector specifying for which subjects to produce the plot. If a character vector, then is
    should contain a subset of the values of the \code{idVar} variable of the \code{newdata} argument of \code{\link{survfitJM}}.}
  \item{fun}{a vectorized function defining a transformation of the survival curve. For example with \code{fun=log} the log-survival curve 
    is drawn.}
  \item{invlink}{a function to transform the fitted values of the longitudinal outcome.}
  \item{conf.int}{logical; if \code{TRUE}, then a pointwise confidence interval is included in the plot.}
  \item{fill.area}{logical; if \code{TRUE} the area defined by the confidence interval of the survival function is put in color.}
  \item{col.area}{the color of the area defined by the confidence interval of the survival function.}
 \item{col.abline,col.points}{the color for the vertical line and the points when \code{include.y} is \code{TRUE}.} 
\item{add.last.time.axis.tick}{logical; if \code{TRUE}, a tick is added in the x-axis for the last available time point for which a
    longitudinal measurement was available.}
  \item{include.y}{logical; if \code{TRUE}, two plots are produced per subject, i.e., the plot of conditional probabilities of survival 
    and a scatterplot of his longitudinal measurements.}
  \item{main}{a character string specifying the title in the plot.}
  \item{xlab}{a character string specifying the x-axis label in the plot.}
  \item{ylab}{a character string specifying the y-axis label in the plot.}
  \item{ylab2}{a character string specifying the y-axis label in the plot, when \code{include.y = TRUE}.}
  \item{lty}{what types of lines to use.}
  \item{col}{which colors to use.}
  \item{lwd}{the thickness of the lines.}
  \item{pch}{the type of points to use.}  
  \item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{par()}.}
  \item{legend}{logical; if \code{TRUE}, a legend is included in the plot.}
  \item{cex.axis.z, cex.lab.z}{the par \code{cex} argument for the axis at side 4, when \code{include.y = TRUE}.}
  \item{xlim}{the par \code{xlim} argument.}
  \item{\dots}{extra graphical parameters passed to \code{plot()}.}
}

\references{
Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for longitudinal and time-to-event data. 
\emph{Biometrics} \bold{67}, 819--829.
}


\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{survfitJM}}}

\examples{
\dontrun{
# we construct the composite event indicator (transplantation or death)
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

# we fit the joint model using splines for the subject-specific 
# longitudinal trajectories and a spline-approximated baseline
# risk function
lmeFit <- lme(log(serBilir) ~ ns(year, 2), data = pbc2,
    random = ~ ns(year, 2) | id)
survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
jointFit <- jointModelBayes(lmeFit, survFit, timeVar = "year")

# we will compute survival probabilities for Subject 2 in a dynamic manner, 
# i.e., after each longitudinal measurement is recorded
ND <- pbc2[pbc2$id == 2, ] # the data of Subject 2
survPreds <- vector("list", nrow(ND))
for (i in 1:nrow(ND)) {
    survPreds[[i]] <- survfitJM(jointFit, newdata = ND[1:i, ])
}

# the default call to the plot method using the first three
# longitudinal measurements
plot(survPreds[[3]])

# we produce the corresponding plot
par(mfrow = c(2, 2), oma = c(0, 2, 0, 2))
for (i in c(1,3,5,7)) {
    plot(survPreds[[i]], estimator = "median", conf.int = TRUE,
        include.y = TRUE, main = paste("Follow-up time:",
            round(survPreds[[i]]$last.time, 1)), ylab = "", ylab2 = "")
}
mtext("log serum bilirubin", side = 2, line = -1, outer = TRUE)
mtext("Survival Probability", side = 4, line = -1, outer = TRUE)
}
}

\keyword{methods}
