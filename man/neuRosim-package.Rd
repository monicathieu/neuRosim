\name{neuRosim-package}
\alias{neuRosim-package}
\alias{neuRosim}
\docType{package}
\title{Functions to generate fMRI data including activated data, noise data and resting state data.
}
\description{The package allows users to generate fMRI time series or 4D data. Some high-level functions are created for fast data generation with only a few arguments and a diversity of functions to define activation and noise. For more advanced users it is possible to use the low-level functions and manipulate the arguments.
}
\details{
\tabular{ll}{
Package: \tab neuRosim\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2010-10-26\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{Marijke Welvaert with contributions from Joke Durnez, Beatrijs Moerkerke, Yves Rosseel and Geert Verdoolaege


Maintainer: Marijke Welvaert <Marijke.Welvaert@Ugent.be>
}
\keyword{ package }
\examples{
## Generate fMRI time series for block design
rm(list=ls())
design <- simprep.act.temp(T=200, ncond=1, onsets=seq(1,200,40), durations=20, TR=2, acc=0.1, hrf="double-gamma")
ts <- sim.fmri.ts(design=design, SNR=1, noise="white")
#plot(ts, type="l")

## Generate fMRI slice for block design with activation in 2 regions
rm(list=ls())
design <- simprep.act.temp(T=200, ncond=1, onsets=seq(1,200,40), durations=20, TR=2, acc=0.1, hrf="double-gamma")
region <- simprep.act.spat(nregio=2, coord=matrix(c(32,15,57,45), ncol=2, byrow=TRUE), ext=c(10,7), form="sphere", fading=TRUE)
out <- sim.fmri.data(design=design, image=region, dim=c(64,64), SNR=1, noise="none")
#plot(out[32,15,], type="l")
}
