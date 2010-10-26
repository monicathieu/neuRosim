\name{canonicalHRF}
\alias{canonicalHRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Double-gamma Haemodynamic reponse function
}
\description{Specifies a double-gamma variate haemodynamic response function for the given time vector and parameters.
}
\usage{
canonicalHRF(x, par = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Time vector in seconds.
}
  \item{par}{Parameter vector of length 5, specifying respectively delay of response (relative to onset), delay of undershoot (relative to onset),
dispersion of response, dispersion of undershoot, scale of undershoot. Default values (i.e. c(6,12,0.9,0.9,0.35)) are used if not specified.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\value{Vector representing the values of the function for the given time vector and parameters.
}
\references{[1] Friston, KJ, Fletcher, P, Josephs, O, Holmes, AP, Rugg, MD and Turner, R (1998). Event-related fMRI: Characterising differential responses. NeuroImage, 7, 30-40.

[2] Glover, GH (1999). Deconvolution of impulse response in event-related BOLD fMRI. NeuroImage, 9, 416-429.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{gammaHRF}}, \code{\link{Balloon}}
}
\examples{
rm(list=ls())
t <- 1:100
out <- canonicalHRF(t, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ activation }% __ONLY ONE__ keyword per line
