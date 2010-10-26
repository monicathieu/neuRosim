\name{gaussnoise.fnc}
\alias{gaussnoise.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Gaussian noise
}
\description{Generates a Gaussian noise dataset with specified dimensions and standard deviation.
}
\usage{
gaussnoise.fnc(dim, sigma, nscan, anatomical = FALSE, template, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{A vector specifying the dimensions of a 2D or 3D array.
}
  \item{sigma}{The standard deviation of the noise.
}
  \item{nscan}{The number of scans in the dataset.
}
  \item{anatomical}{Logical indicating if an anatomical structure should be incorporated.
}
  \item{template}{If anatomical==TRUE, an array representing the anatomical structure with dimensions equal to dim.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\details{The function generates random Gaussian noise. The result is a noise
 array with specified dimensions and desired standard deviation.
}
\value{An array containing the noise with dimensions specified in dim and nscan.
}

\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AR1noise.fnc}}, \code{\link{lowfreq.fnc}}, \code{\link{physnoise.fnc}}, \code{\link{tasknoise.fnc}}, \code{\link{spatnoise.fnc}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
sigma <- 5
nscan <- 100
out <- gaussnoise.fnc(d, sigma, nscan, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ noise }% __ONLY ONE__ keyword per line
