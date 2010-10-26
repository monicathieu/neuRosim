\name{gammaHRF}
\alias{gammaHRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single Gamma Haemodynamic response function.
}
\description{Specifies a Gamma variate haemodynamic response function for the given time vector and FWHM.
}
\usage{
gammaHRF(x, FWHM = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Time vector in seconds.
}
  \item{FWHM}{Full Width Half Maximum of the Gamma variate function. Default value (i.e. 4) is used if not specified.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\value{Vector representing the values of the function for the given time vector and FWHM.
}
\references{Buxton, RB, Uludag, K, Dubowitz, DJ and Liu, TT (2004). Modeling the hemodynamic response to brain activation. NeuroImage, 23, S220-S233.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{canonicalHRF}}, \code{\link{Balloon}}
}
\examples{
rm(list=ls())
t <- 1:100
out <- gammaHRF(t, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level}
\keyword{ activation }% __ONLY ONE__ keyword per line
