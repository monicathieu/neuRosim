\name{simprep.act.temp}
\alias{simprep.act.temp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare temporal structure of the data
}
\description{Prepare a list defining the necessary parameters to specify the temporal structure of the activation data.
}
\usage{
simprep.act.temp(T, ncond, nregio = NULL, onsets, durations, TR, acc, hrf = c("gamma", "double-gamma", "Balloon"), par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{Duration of the experiment.
}
  \item{ncond}{Number of conditions (i.e. columns in the design matrix).
}
  \item{nregio}{Number of regions. If not specified, it is assumed that all regions have the same design matrix.
}
  \item{onsets}{Vector, matrix or array representing the onsets of the stimulus in seconds.
}
  \item{durations}{Vector, matrix or array representing the durations of the stimulus in seconds.
}
  \item{TR}{Repetition time in seconds.
}
  \item{acc}{Microtime resolution in seconds.
}
  \item{hrf}{Haemodynamic response function (double-gamma is default)
}
  \item{par}{Vector, matrix or array representing the parameters of the haemodynamic response function.
}
}
\details{If each region has its own design matrix, onsets and durations should be an array with the second dimension corresponding to the number of conditions and the third dimension corresponding to the number of regions.
}
\value{A list with the necessary arguments to be used in \code{\link{sim.fmri.data}} or \code{\link{sim.fmri.ts}}.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sim.fmri.data}}, \code{\link{sim.fmri.ts}}, \code{\link{simprep.act.spat}}, \code{\link{specify.region}}
}
\examples{
rm(list=ls())
ncond <- 2
os <- matrix(c(20,60,15,35), ncol=ncond)
d <- c(20, 10)
T <- 80
TR <- 2
acc <- 0.1
out <- simprep.act.temp(T, ncond, onsets=os, durations=d, TR=TR, acc=acc, hrf="Balloon")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ high-level }
