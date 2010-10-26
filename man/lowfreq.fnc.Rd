\name{lowfreq.fnc}
\alias{lowfreq.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate low frequency drift
}
\description{Generates a low-frequency drift dataset with specified dimensions and frequency.
}
\usage{
lowfreq.fnc(dim, freq = NULL, nscan, TR, anatomical = FALSE, template, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{A vector specifying the dimensions of a 2D or 3D array.
}
  \item{freq}{The frequency of the drift in seconds. Default value (i.e. 128) is used when not specified.
}
  \item{nscan}{The number of scans in the dataset.
}
  \item{TR}{The repetition time in seconds.
}
  \item{anatomical}{Logical indicating if an anatomical structure should be incorporated.
}
  \item{template}{If anatomical==TRUE, an array representing the anatomical structure with dimensions equal to dim.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\details{The function generates low-frequency drift based on a basis set of cosine functions. The result is an array with specified dimensions
and frequency.
}
\value{An array containing the drift with dimensions specified in dim.
}
\references{Friston et al. (2007). Statistical Parametric Mapping: The analysis of functional brain images. Academic Press.
}
\author{Y. Rosseel, M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AR1noise.fnc}}, \code{\link{gaussnoise.fnc}}, \code{\link{physnoise.fnc}}, \code{\link{tasknoise.fnc}}, \code{\link{spatnoise.fnc}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
freq <- 80
nscan <- 100
TR <- 2
out <- lowfreq.fnc(d, freq, nscan, TR, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ noise}% __ONLY ONE__ keyword per line
