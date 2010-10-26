\name{sim.restingstate.ts}
\alias{sim.restingstate.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate fMRI resting state time series
}
\description{Synthesizes a single time series x representing resting state activity. The fluctuation frequencies f are limited to 
a square passband 0.01 Hz <= f <= 0.1 Hz. TR is the repetition time (needed to compute the passband limits), expressed in seconds. 
N is the required number of samples (needs not be a power of 2).
}
\usage{
sim.restingstate.ts(nscan, TR, SNR, noise = c("none", "white", "temporal", "low-frequency", "physiological", "task-related", "mixture"), temp = c("AR(1)", "state-space"), weights, verbose = TRUE, rho = NULL, lowfreq = NULL, heartfreq = NULL, respfreq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nscan}{Number of scans.
}
  \item{TR}{Repetition time in seconds.
}
  \item{SNR}{Signal-to-noise ratio of the time series.
}
  \item{noise}{Type of noise (white is default).
}
  \item{temp}{If noise=temporal, type of temporal noise (AR(1) is default).
}
  \item{weights}{If noise=mixture, vector of weights to specify the fraction of the noise components.
}
  \item{verbose}{Logical indicating if warnings should be returned.
}
  \item{rho}{If noise=temporal or noise=mixture, value of autocorrelation coefficient.
}
  \item{lowfreq}{If noise=low-frequency or noise=mixture, frequency of the low-frequency drift in seconds.
}
  \item{heartfreq}{If noise=physiological or noise=mixture, frequency of heart rate in Hz.
}
  \item{respfreq}{If noise=physiological or noise=mixture, frequency of respiratory rate in Hz.
}
}
\value{A vector representing the resting state time series
}
\references{[1] C.G. Fox, Computers & Geoscience, Vol. 13, pp. 369-374, 1987.
 [2] M. Fukunaga, Magnetic Resonance Imaging, Vol. 24, pp. 979-992, 2006.
}
\author{J. Durnez, G. Verdoolaege, M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sim.fmri.ts}}
}
\examples{
rm(list=ls())
out <- sim.restingstate.ts(nscan=50, TR=2, SNR=1, noise="none")
#plot(out, type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ high-level }
