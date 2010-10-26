\name{specify.region}
\alias{specify.region}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate activation image
}
\description{Generates an image with activated regions for specified dimensions. The regions are defined by there center and radius or can be entered manually.
}
\usage{
specify.region(dim, coord, ext = NULL, form = c("cube", "sphere", "manual"), fading = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{Dimensions of the image space.
}
  \item{coord}{Coordinates of the activated region, if form=cube or sphere, the coordinates represent the center of the region, if form=manual, the coordinates should be in matrix form, where the rows represent the voxels and the columns the x-y-z coordinates.
}
  \item{ext}{If form=cube or sphere, the distance in voxels from the center of the region to the edge.
}
  \item{form}{The form of the activated region. Default is "cube".
}
  \item{fading}{Logical indicating if the center should be more activated than the edges.
}
}
\value{An array representing the activation image with specified regions.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{specify.region}},\code{\link{gammaHRF}},\code{\link{canonicalHRF}},\code{\link{Balloon}}
}
\examples{
rm(list=ls())
d <- c(10,10,10)
coord <- c(3,3,3)
ext <- 1
out <- specify.region(d, coord, ext, form="sphere", fading=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
