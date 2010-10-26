\name{physnoise.fnc}
\alias{physnoise.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate physiological noise
}
\description{Generates a physiological noise dataset with specified dimensions and standard deviation. The physiological noise is defined as noise caused by heart beat and respiratory rate.
}
\usage{
physnoise.fnc(dim, sigma, nscan, TR, freq.heart = NULL, freq.resp = NULL, anatomical = FALSE, template, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{A vector specifying the dimensions of a 2D or 3D array.
}
  \item{sigma}{The standard deviation of the noise.
}
  \item{nscan}{The number of scans in the dataset.
}
  \item{TR}{The repetition time in seconds.
}
  \item{freq.heart}{The frequency in Hz of the heart beat. Default value (i.e. 1.17) is used if not specified.
}
  \item{freq.resp}{The frequency in Hz of the respiratory rate. Default value (i.e. 0.2) is used if not specified.
}
  \item{anatomical}{Logical indicating if an anatomical structure should be incorporated.
}
  \item{template}{If anatomical==TRUE, an array representing the anatomical structure with dimensions equal to dim.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\details{The function generates physiological noise. Heart beat and respiratory rate are defined as sine and cosine functions with specified frequencies. Additional Gaussian noise creates variability over voxels. The result is a noise dataset with specified dimensions and desired standard deviation.
}
\value{An array containing the noise with dimensions specified in dim and nscan.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AR1noise.fnc}}, \code{\link{lowfreq.fnc}}, \code{\link{gaussnoise.fnc}}, \code{\link{tasknoise.fnc}}, \code{\link{spatnoise.fnc}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
sigma <- 5
nscan <- 100
TR <- 2
out <- physnoise.fnc(d, sigma, nscan, TR, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ noise }% __ONLY ONE__ keyword per line
