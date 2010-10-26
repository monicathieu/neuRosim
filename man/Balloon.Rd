\name{Balloon}
\alias{Balloon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Balloon model
}
\description{Generates the BOLD signal based on the Balloon model for a given stimulus vector and microtime resolution.
}
\usage{
Balloon(stim, T, acc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stim}{Vector representing the presence/absence (1-0 coding) of a stimulus/activation in seconds.
}
  \item{T}{Total duration of stimulus vector in seconds.
}
  \item{acc}{Microtime resolution of stimulus vector in seconds.
}
}
\value{Vector representing the values of the BOLD signal for the given stimulus vector.
}
\references{Buxton, RB, Uludag, K, Dubowitz, DJ and Liu, TT (2004). Modeling the hemodynamic response to brain activation. NeuroImage, 23, S220-S233.
}
\author{G. Verdoolaege, M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{This function is still in developing fase. Additional functionalities to manipulate the parameters of the Balloon model are not yet implemented.}
\seealso{\code{\link{canonicalHRF}}, \code{\link{gammaHRF}}
}
\examples{
# NOT RUN (takes a couple of seconds)
#s <- rep(rep(0,10), rep(1,10), 5)
#T <- 100
#it <- 0.1
#out <- Balloon(s, T, it)
# END(NOT RUN)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ activation}% __ONLY ONE__ keyword per line
