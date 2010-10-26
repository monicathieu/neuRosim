\name{tasknoise.fnc}
\alias{tasknoise.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate task-related noise
}
\description{Generates a Gaussian noise dataset with specified dimensions and standard deviation only when a task is performed or activation is present.
}
\usage{
tasknoise.fnc(dim, act.image, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{A vector specifying the dimensions of a 2D or 3D array.
}
  \item{act.image}{Array defining where and when activation is present.
}
  \item{sigma}{Standard deviation of the noise.
}
}
\details{The function generates random Gaussian noise for those voxels in the dataset that show activation. The result is a noise
 array with specified dimensions and desired standard deviation.
}
\value{An array containing the noise with dimensions specified in dim and nscan.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AR1noise.fnc}}, \code{\link{lowfreq.fnc}}, \code{\link{physnoise.fnc}}, \code{\link{gaussnoise.fnc}}, \code{\link{spatnoise.fnc}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
sigma <- 5
nscan <- 100
act <- array(rep(0, prod(d)*nscan), dim=c(d,nscan))
act[2:4,2:4,2:4,c(20:30,40:50,60:70)] <- 1
out <- tasknoise.fnc(d, act, sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ noise }% __ONLY ONE__ keyword per line
