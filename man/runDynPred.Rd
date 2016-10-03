\name{runDynPred}

\alias{runDynPred}



\title{ Shiny Application for Dynamic Predictions }

\description{
    This function loads the shiny package and runs the application for calculating dynamic predictions using package JMbayes.
}

\usage{
runDynPred(type = c("JM", "lme"))
} 

\arguments{
    \item{type}{character string indicating whether dynamic predictions are based on 
        joint models or mixed models alone.}
}

\value{
No value returned.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
\dontrun{
runDynPred()
}
}

\keyword{misc}

