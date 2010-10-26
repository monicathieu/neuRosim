\name{data.read.nifti}
\alias{data.read.nifti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read data from NIFTI files
}
\description{This function reads the data from a nifti file.
}
\usage{
data.read.nifti(filename = "", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{The name of the nifti file.
}
  \item{verbose}{Logical indicating if processing information should be displayed.
}
}
\value{An array representing the voxel values.
}
\references{http://nifti.nimh.nih.gov, http://niftilib.sourceforge.net
}
\author{Y. Rosseel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{data.write.nifti}}
}
\examples{
## Not run:
#data.read.nifti("file", verbose=FALSE)
## End(Not run)
}
\keyword{ low-level }
