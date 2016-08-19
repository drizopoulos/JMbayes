\name{IndvPred_lme}

\alias{IndvPred_lme}

\title{Individualized Predictions from Linear Mixed Models}

\description{
  Calculates subject-specific predictions for new subjects from a linear mixed model.
}

\usage{
IndvPred_lme(lmeObject, newdata, timeVar, times = NULL, M = 200L,
    interval = c("confidence", "prediction"), level = 0.95, 
    return_data = FALSE, seed = 1L)
}

\arguments{
  \item{lmeObject}{an object inheriting from class \code{lme}.}
  \item{newdata}{a data frame in which to look for variables with which to predict.}
  \item{timeVar}{a character string specifying the time variable in the linear mixed 
        model.}
  \item{interval}{a character string indicating what type of intervals should be computed.}
  \item{level}{a numeric scalar denoting the tolerance/confidence level.}
  \item{times}{a numeric vector denoting the time points for which we wish to compute the 
        subject-specific predictions after the last available measurement provided in 
        \code{newdata}. Default is a sequence of 100 equally spaced time points 
        from the smallest to the largest follow-up time of all subjects.}
  \item{M}{numeric scalar denoting the number of Monte Carlo samples. 
        See \bold{Details}.}
  \item{return_data}{logical; if \code{TRUE} the data frame supplied in
        \code{newdata} is returned augmented with the outputs of the function.}
  \item{seed}{numeric scalar, the random seed used to produce the results.}
}

\value{
  If \code{return_data = TRUE}, a the data frame \code{newdata} with extra rows for the
  time points at which predictions were calculated, and extra columns with the predictions
  and the limits of the pointwise confidence intervals.
  
  If \code{return_data = FALSE}, a list with components
  \item{times_to_pred}{time points at which predictions were calculated.}
  \item{predicted_y}{the predictions.}
  \item{low}{the lower limits of the pointwise confidence intervals.}
  \item{upp}{the upper limits of the pointwise confidence intervals.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{predict.JMbayes}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * ns(year, 2), data = pbc2, 
        random = ~ ns(year, 2) | id)
        
DF <- IndvPred_lme(fitLME, newdata = subset(pbc2, id %in% c(2, 25)), timeVar = "year", 
    M = 500, return_data = TRUE)
    
require(lattice)
xyplot(pred + low + upp ~ year | id, data = DF,
    type = "l", col = c(2,1,1), lty = c(1,2,2), lwd = 2,
    ylab = "Average log serum Bilirubin")
}
}
\keyword{methods}
