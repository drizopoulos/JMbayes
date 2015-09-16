\name{gt}

\alias{dgt}
\alias{pgt}
\alias{qgt}
\alias{rgt}



\title{ The Generalized Student's t Distribution }

\description{
    Density, distribution function, quantile function and random generation for the generalized Student's t distribution.
}

\usage{
dgt(x, mu = 0, sigma = 1, df = stop("no df arg"), log = FALSE)

pgt(q, mu = 0, sigma = 1, df = stop("no df arg")) 

qgt(p, mu = 0, sigma = 1, df = stop("no df arg"))

rgt(n, mu = 0, sigma = 1, df = stop("no df arg"))
} 

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{a numeric scalar denoting the number of observations.}
  \item{mu}{a vector of means.}
  \item{sigma}{a vector of standard deviations.}
  \item{log}{logical; if \code{TRUE} the density is computed in the log scale.}
  \item{df}{a numeric scalar denoting the degrees of freedom.}
}


\value{
 \code{dgt} gives the density, \code{pgt} gives the distribution function, 
 \code{qgt} gives the quantile function, and \code{rgt} generates random deviates.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
x <- rnorm(10, mean = 10, sd = 3)
dgt(x, mu = 10, sigma = 3, df = 4)
rgt(10, mu = 10, sigma = 3, df = 4)
}

\keyword{misc}

