\name{spatnoise.fnc}
\alias{spatnoise.fnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate spatially correlated noise
}
\description{Generates a spatially correlated noise dataset with specified dimensions and standard deviation.
}
\usage{
spatnoise.fnc(dim, sigma, nscan, method = c("corr", "gammaRF", "gausRF"), rho = NULL, FWHM = NULL, gamma.shape = NULL, gamma.rate = NULL, anatomical = FALSE, template, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{A vector specifying the dimensions of a 2D or 3D array.
}
  \item{sigma}{The standard deviation of the noise.
}
  \item{nscan}{The number of scans in the dataset.
}
  \item{method}{Method specifying the type of spatial correlation. Default is "corr".
}
  \item{rho}{If method=="corr", the value of the autocorrelation coefficient. Default value (i.e. 0.75) is used if not specified.
}
  \item{FWHM}{If method=="gammaRF" or "gausRF", the full-width-half-maximum of the Gaussian kernel. Default value (i.e. 4) is used if not specified.
}
  \item{gamma.shape}{If method=="gammaRF", the shape parameter of the Gamma distribution. Default value (i.e. 6) is used if not specified.
}
  \item{gamma.rate}{If method=="gammaRF", the shape parameter of the Gamma distribution. Default value (i.e. 1) is used if not specified.
}
  \item{anatomical}{Logical indicating if an anatomical structure should be incorporated.
}
  \item{template}{If anatomical==TRUE, an array representing the anatomical structure with dimensions equal to dim.
}
  \item{verbose}{Logical indicating if warnings should be printed.
}
}
\details{The function generates spatially correlated noise. When method=="corr", AR(1) voxelwise correlations are introduced. 
If method=="gausRF" of "gammaRF", respectively a Gaussian Random Field or a Gamma Random Field is created. The result is a noise array with specified dimensions and desired standard deviation.
}
\value{An array containing the noise with dimensions specified in dim and nscan.
}
\author{J. Durnez, B. Moerkerke, M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AR1noise.fnc}}, \code{\link{lowfreq.fnc}}, \code{\link{physnoise.fnc}}, \code{\link{tasknoise.fnc}}, \code{\link{gaussnoise.fnc}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
sigma <- 5
nscan <- 100
rhospat <- 0.7
out <- spatnoise.fnc(d, sigma, nscan, method="corr", rho=rhospat, verbose=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
\keyword{ noise }% __ONLY ONE__ keyword per line
