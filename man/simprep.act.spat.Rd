\name{simprep.act.spat}
\alias{simprep.act.spat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare spatial structure of the data
}
\description{Prepare a list defining the necessary parameters to specify the spatial structure of the activation data.
}
\usage{
simprep.act.spat(nregio, coord, ext = NULL, form = c("cube", "sphere", "manual"), fading = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nregio}{Number of activated regions.
}
  \item{coord}{Matrix of coordinates, the columns represent the xyz-coordinates.
}
  \item{ext}{If form=cube or sphere, the distance between the center and the edge, if form=manual, the number of voxels in each region.
}
  \item{form}{The form of the activated regions.
}
  \item{fading}{Logical indicating if the center of the region should be more activated than the edges.
}
}
\value{A list with the necessary arguments to be used in \code{\link{sim.fmri.data}}.
}
\author{M. Welvaert
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sim.fmri.data}}, \code{\link{simprep.act.temp}}, \code{\link{specify.region}}
}
\examples{
rm(list=ls())
coord <- matrix(c(3,3,3,6,6,6), ncol=3, byrow=TRUE)
ext <- c(1,2)
out <- simprep.act.spat(2, coord, ext, form="cube", fading=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ high-level }
