\name{stim.boxcar}
\alias{stim.boxcar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a stimulus boxcar function.
}
\description{Generates a stimulus boxcar vector for the specified time duration and microtime resolution based on the user-defined onsets and durations.
}
\usage{
stim.boxcar(T, onset, duration, acc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Total time of the design in seconds.
}
  \item{onset}{Vector representing the onsets of the stimulus in seconds.
}
  \item{duration}{Vector representing the durations of the stimulus in seconds.
}
  \item{acc}{Microtime resolution in seconds.
}
}
\details{If duration is a single number, it is assumed that all stimulus onsets have the same duration. 
}
\value{A vector in microtime resolution specifying the stimulus boxcar function in 1-0 coding.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{specify.design}}
}
\examples{
rm(list=ls())
T <- 100
os <- c(1, 21, 41, 61, 81)
d <- 10
out <- stim.boxcar(T, os, d, 0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
