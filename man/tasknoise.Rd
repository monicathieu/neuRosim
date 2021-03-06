% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noise.R
\name{tasknoise}
\alias{tasknoise}
\title{Generate task-related noise}
\usage{
tasknoise(act.image, sigma, type = c("gaussian", "rician"), vee = 1,
  verbose = T)
}
\arguments{
\item{act.image}{Array defining where and when activation is present.}

\item{sigma}{Standard deviation of the noise.}

\item{type}{Distribution of task-related noise. Default is gaussian.}

\item{vee}{If \code{type=="rician"}, the non-centrality parameter of the distribution.}

\item{verbose}{Logical indicating if warnings should be printed.}
}
\value{
An array containing the noise.
}
\description{
Generates a Gaussian noise dataset with specified dimensions and standard deviation
only when a task is performed or activation is present.
}
\details{
The function generates random Gaussian noise for those voxels in the dataset
that show activation. The result is a noise array with specified dimensions and
desired standard deviation.
}
\examples{
d <- c(10,10,10)
nscan <- 100
act <- array(rep(0, prod(d)*nscan), dim=c(d,nscan))
act[2:4,2:4,2:4,c(20:30,40:50,60:70)] <- 1
out <- tasknoise(act, sigma = 5)
\dontshow{rm(d,nscan,act,out)}

}
\seealso{
\code{\link{temporalnoise}}, \code{\link{lowfreqdrift}}, \code{\link{physnoise}},
\code{\link{systemnoise}}, \code{\link{spatialnoise}}
}
\keyword{low-level}
\keyword{noise}
