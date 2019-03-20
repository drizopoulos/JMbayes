\name{tve}

\alias{tve}

\title{ Time-Varying Effects using P-splines }

\description{
    A B-spline expansion of the input variables to be used for a time-varying effect 
    in the specification of joint model.
}

\usage{
tve(x, df = NULL, knots = NULL, ord = 3)
} 

\arguments{
  \item{x}{a numeric input variable.}
  \item{df}{integer denoting the degrees of freedom.}
  \item{knots}{a numeric vector of knots.}
  \item{ord}{an integer denoting the order of the spline.}
}


\value{
  an object of class \code{tve}.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\keyword{regression}

