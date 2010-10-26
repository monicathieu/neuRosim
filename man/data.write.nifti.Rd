\name{data.write.nifti}
\alias{data.write.nifti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write data to NIFTI files.
}
\description{This function write data arrays to a NIFTI file.
}
\usage{
data.write.nifti(V, output.file, check = 1, datatype = 16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{Data array
}
  \item{output.file}{Name of the output file.
}
  \item{check}{Check for previous existence of filename (existence is an error condition).
}
  \item{datatype}{Bit size of the data points.
}
}
\references{http://nifti.nimh.nih.gov, http://niftilib.sourceforge.net
}
\author{Y. Rosseel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{data.read.nifti}}
}
\examples{
## Not run:
#data <- array(rnorm(64*64*10), dim=c(64,64,1,10))
#data.write.nifti(data, "filename")
## End(Not run)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ low-level }
