\name{specify.design}
\alias{specify.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate design matrix.
}
\description{Generates a design matrix to be used as a model for the simulated activation.
}
\usage{
specify.design(ncond, onsets, durations, T, TR, acc, conv = c("none", "gamma", "double-gamma", "Balloon"), cond.names = NULL, par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncond}{Number of conditions, i.e. columns in the design matrix.
}
  \item{onsets}{Vector or matrix representing the onsets in seconds, ncol should match ncond.
}
  \item{durations}{Vector or matrix representing the durations in seconds, ncol should match ncond.
}
  \item{T}{Duration of the experiment in seconds.
}
  \item{TR}{Repetition time in seconds.
}
  \item{acc}{Microtime resolution in seconds.
}
  \item{conv}{Should the design matrix be convoluted, default is none. 
}
  \item{cond.names}{Optional names for the conditions.
}
  \item{par}{Parameters of the haemodynamic response function. See \code{\link{gammaHRF}} and \code{\link{canonicalHRF}} for more details.
}
}
\value{A matrix specifying the design.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{specify.region}},\code{\link{gammaHRF}},\code{\link{canonicalHRF}},\code{\link{Balloon}}
}
\examples{
rm(list=ls())
ncond <- 2
os <- matrix(c(20,60,15,35), ncol=ncond)
d <- c(20, 10)
T <- 80
TR <- 2
acc <- 0.1
out <- specify.design(ncond, os, d, T, TR, acc, conv="double-gamma")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level}
