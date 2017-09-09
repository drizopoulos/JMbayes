\name{DerivSplines}

\alias{dns}
\alias{dbs}
\alias{ins}
\alias{ibs}



\title{ Derivatives and Integrals of B-splines and Natural Cubic splines }

\description{
    Numerical derivatives and integrals of functions \code{bs()} and \code{ns()} at their first argument.
}

\usage{
dns(x, df = NULL, knots = NULL, intercept = FALSE, 
    Boundary.knots = range(x), eps = 1e-03)

dbs(x, df = NULL, knots = NULL, intercept = FALSE, 
    Boundary.knots = range(x), eps = 1e-03)

ins(x, df = NULL, knots = NULL, intercept = FALSE, 
    Boundary.knots = range(x), from = 0, weight.fun = NULL, 
    integrand.fun = NULL, \dots)

ibs(x, df = NULL, knots = NULL, intercept = FALSE, 
    Boundary.knots = range(x), from = 0, weight.fun = NULL, 
    integrand.fun = NULL, \dots)
} 

\arguments{
  \item{x, df, knots, intercept, Boundary.knots}{see the help pages of functions \code{ns()} and \code{bs()}.}
  \item{eps}{a numeric scalar denoting the step length for the central difference approximation, which
            calculates the derivative.}
  \item{from}{a numeric scalar denoting the lower limit of the integral.}
  \item{weight.fun}{a function to be applied as weights.}
  \item{integrand.fun}{a function to be applied in the integrand.}
  \item{\dots}{extra arguments passed to \code{weight.fun}.}
}


\value{
  an object of class \code{dns}, \code{dbs}, \code{ins} or \code{ibs}.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
x <- rnorm(10)
dns(x, df = 4)
ins(x, df = 4)
}

\keyword{regression}

