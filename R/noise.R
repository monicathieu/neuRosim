#' Generate system noise
#' 
#' Generates a system noise dataset with specified dimensions and standard deviation.
#' The noise can be either Gaussian or Rician distributed.
#' 
#' @export
#' 
#' @param dim A vector specifying the dimensions of the image.
#' @param nscan The number of images in the dataset.
#' @param type Distribution of system noise. Default is gaussian.
#' @param sigma The standard deviation of the noise.
#' @param vee If \code{type=="rician"}, the non-centrality parameter of the distribution .
#' @param template An array representing the anatomical structure or mask with dimensions equal to dim.
#' @param verbose Logical indicating if warnings should be printed.
#' @return An array containing the noise with dimensions specified in \code{dim} and \code{nscan}.
#' 
#' @seealso \code{\link{temporalnoise}}, \code{\link{lowfreqdrift}}, \code{\link{physnoise}},
#' \code{\link{tasknoise}}, \code{\link{spatialnoise}}
#' 
#' @examples
#' d <- c(10,10,10)
#' sigma <- 5
#' nscan <- 100
#' out <- systemnoise(d, nscan, type="rician", sigma, verbose=FALSE)
#' \dontshow{rm(d,sigma,nscan,out)}
#' 
#' @keywords low-level noise

systemnoise <- function(dim, nscan, type=c("gaussian","rician"), sigma, vee=1, template, verbose=TRUE){
  
  if(length(dim)>3) stop("Image space with more than three dimensions is not supported.")
  
  if(missing(type))	type="gaussian"
  
  if (verbose) cat("System noise type: ", type, "\n")
  
  if(type=="gaussian"){
    noise <- array(rnorm(prod(dim)*nscan, 0, sigma), dim=c(dim,nscan))
  }else{ 
    if(type=="rician"){
      noise <- array(rrice(prod(dim)*nscan, vee=vee, sigma=sigma), dim=c(dim,nscan))
    } else {
      stop("Specified type of system noise is unknown. Type should be gaussian or rician.")
    }
  }
  
  if(!missing(template)){
    if(length(dim(template))>3){
      stop("Template should be a 2D or 3D array.")
    }
    template.time <- array(rep(template,nscan), dim=c(dim,nscan))
    ix <- which(template.time!=0)
    noise[-ix] <- 0
  }
  
  return(noise)
}

#' Generate temporally correlated noise
#' 
#' Generates an autoregressive noise dataset with specified dimensions and standard deviation.
#' 
#' @export
#' @importFrom stats arima.sim
#' 
#' @param dim A vector specifying the dimensions of a 2D or 3D array.
#' @param nscan The number of scans in the dataset.
#' @param sigma The standard deviation of the noise.
#' @param rho The autocorrelation coefficients. The length of the vector determines the order of the autoregressive model.
#' @param template An array representing the anatomical structure or mask with dimensions equal to dim.
#' @param verbose Logical indicating if warnings should be printed.
#' @return An array containing the noise with dimensions specified in \code{dim}.
#' 
#' @seealso \code{\link{systemnoise}}, \code{\link{lowfreqdrift}}, \code{\link{physnoise}},
#' \code{\link{tasknoise}}, \code{\link{spatialnoise}}
#' 
#' @examples
#' d <- c(10,10,10)
#' sigma <- 5
#' nscan <- 100
#' rho <- c(0.3,-0.7)
#' out <- temporalnoise(d, nscan, sigma, rho, verbose=FALSE)
#' \dontshow{rm(d,sigma,nscan,rho,out)}
#' 
#' @keywords low-level noise

temporalnoise <- function(dim, nscan, sigma, rho=0.2, template, verbose=TRUE){
  
  if(length(dim)>3){
    stop("Image space with more than three dimensions is not supported.")
  }
  
  if (verbose) cat("Temporal noise: using rho = ", rho, "\n")
  
  array.dim <- c(dim,rep(1,3-length(dim)),nscan)
  noise.array <- array(0, dim=array.dim)
  
  for(i in 1:array.dim[1]){
    for(j in 1:array.dim[2]){
      for(k in 1:array.dim[3]){
        noise.array[i,j,k,] <- arima.sim(list(ar=rho, ma=0), sd=sigma, n=nscan)
      }
    }
  }
  
  noise <- array(c(noise.array), dim=c(dim, nscan))
  
  if(!missing(template)){
    if(length(dim(template))>3){
      stop("Template should be a 2D or 3D array.")
    }
    template.time <- array(rep(template,nscan), dim=c(dim,nscan))
    ix <- which(template.time!=0)
    noise[-ix] <- 0
  }
  
  return(noise)
}

#' Generate low frequency drift
#' 
#' Generates a low-frequency drift dataset with specified dimensions and frequency.
#' 
#' @export
#' 
#' @param dim A vector specifying the dimensions of the image.
#' @param freq The frequency of the drift in seconds.
#' @param nscan The number of scans in the dataset.
#' @param TR The repetition time in seconds.
#' @param template An array representing the anatomical structure or mask with dimensions equal to \code{dim}.
#' @param verbose Logical indicating if warnings should be printed.
#' @return An array containing the drift with dimensions specified in \code{dim}.
#' 
#' @details The function generates low-frequency drift based on a basis set of cosine functions.
#' The result is an array with specified dimensions and frequency.
#' 
#' @references Friston et al. (2007). Statistical Parametric Mapping: The analysis of functional brain images. Academic Press.
#' 
#' @seealso \code{\link{temporalnoise}}, \code{\link{systemnoise}}, \code{\link{physnoise}},
#' \code{\link{tasknoise}}, \code{\link{spatialnoise}}
#' 
#' @examples
#' d <- c(10,10,10)
#' freq <- 80
#' nscan <- 100
#' TR <- 2
#' out <- lowfreqdrift(d, freq, nscan, TR, verbose=FALSE)
#' \dontshow{rm(d,freq,nscan,TR,out)}
#' 
#' @keywords low-level noise

lowfreqdrift <- function(dim, freq=128, nscan, TR, template, verbose=TRUE){
  
  if (verbose) cat("Using low-frequency drift period = ", freq, "\n")
  
  n <- floor( 2*(nscan*TR)/freq + 1 )
  if(n<3){
    stop("Number of basis functions is too low. Lower frequency or longer scanning time should be provided.")
  }
  drift.base <- spm_drift(nscan, n)[,-1]  
  drift.image <- array(rep(1, prod(dim)), dim=c(dim))
  drift.array <- drift.image %o% rowSums(drift.base)
  
  if(!missing(template)){
    if(length(dim(template))>3){
      stop("Template should be a 2D or 3D array.")
    }
    template.time <- array(rep(template,nscan), dim=c(dim,nscan))
    ix <- which(template.time!=0)
    drift.array[-ix] <- 0
  }
  
  return(drift.array)
}

spm_drift <- function(N, K) {
  # Helper function for lowfreqdrift
  n <- 0:(N-1)
  C <- matrix(0, nrow=N, ncol=K)
  
  C[,1] <- 1/sqrt(N)
  for(k in 2:K) {
    C[,k] = sqrt(2/N)*10*cos(pi*(2*n+1)*(k-1)/(2*N))
  }
  return(C)
}

#' Generate physiological noise
#' 
#' Generates a physiological noise dataset with specified dimensions and standard deviation.
#' The physiological noise is defined as noise caused by heart beat and respiratory rate.
#' 
#' @export
#' @importFrom stats sd rnorm
#' 
#' @param dim A vector specifying the dimensions of the image.
#' @param nscan The number of scans in the dataset.
#' @param TR The repetition time in seconds.
#' @param sigma The standard deviation of the noise.
#' @param freq.heart The frequency in Hz of the heart beat.
#' @param freq.resp The frequency in Hz of the respiratory rate.
#' @param template An array representing the anatomical structure or mask with dimensions equal to dim.
#' @param verbose Logical indicating if warnings should be printed.
#' @return An array containing the noise with dimensions specified in \code{dim} and \code{nscan}.
#' 
#' @details The function generates physiological noise. Heart beat and respiratory rate are defined
#' as sine and cosine functions with specified frequencies. Additional Gaussian noise creates
#' variability over voxels. The result is a noise dataset with specified dimensions and desired
#' standard deviation.
#' 
#' @seealso \code{\link{temporalnoise}}, \code{\link{lowfreqdrift}}, \code{\link{systemnoise}},
#' \code{\link{tasknoise}}, \code{\link{spatialnoise}}
#' 
#' @examples
#' d <- c(10,10,10)
#' sigma <- 5
#' nscan <- 100
#' TR <- 2
#' out <- physnoise(d, nscan, TR, sigma, verbose=FALSE)
#' \dontshow{rm(d,sigma,nscan,TR,out)}
#' 
#' @keywords low-level noise

physnoise <- function(dim, nscan, TR, sigma, freq.heart=1.17, freq.resp=0.2, template, verbose=TRUE){
  
  if (verbose) cat("Using cardiac frequency = ", freq.heart," Hz\n",
                   "Using respiratory frequency = ", freq.resp, " Hz\n")
  
  if(length(dim)>3){
    stop("Image space with more than three dimensions is not supported.")
  }
  
  HB <- 2*pi*freq.heart*TR
  RR <- 2*pi*freq.resp*TR
  t <- 1:nscan
  
  HRdrift <- sin(HB*t) + cos(RR*t)
  sigma.HR <- sd(HRdrift)
  HRweight <- sigma/sigma.HR
  
  noise <- array(rnorm(prod(dim)*nscan, 0, 1), dim=c(dim, nscan)) + HRweight*HRdrift
  
  if(!missing(template)){
    if(length(dim(template))>3){
      stop("Template should be a 2D or 3D array.")
    }
    template.time <- array(rep(template,nscan), dim=c(dim,nscan))
    ix <- which(template.time!=0)
    noise[-ix] <- 0
  }
  
  return(noise)
}

#' Generate task-related noise
#' 
#' Generates a Gaussian noise dataset with specified dimensions and standard deviation
#' only when a task is performed or activation is present.
#' 
#' @export
#' @importFrom stats rnorm sd
#' 
#' @param act.image Array defining where and when activation is present.
#' @param sigma Standard deviation of the noise.
#' @param type Distribution of task-related noise. Default is gaussian.
#' @param vee If \code{type=="rician"}, the non-centrality parameter of the distribution.
#' @param verbose Logical indicating if warnings should be printed.
#' @return An array containing the noise.
#' 
#' @details The function generates random Gaussian noise for those voxels in the dataset
#' that show activation. The result is a noise array with specified dimensions and
#' desired standard deviation.
#' 
#' @seealso \code{\link{temporalnoise}}, \code{\link{lowfreqdrift}}, \code{\link{physnoise}},
#' \code{\link{systemnoise}}, \code{\link{spatialnoise}}
#' 
#' @examples
#' d <- c(10,10,10)
#' nscan <- 100
#' act <- array(rep(0, prod(d)*nscan), dim=c(d,nscan))
#' act[2:4,2:4,2:4,c(20:30,40:50,60:70)] <- 1
#' out <- tasknoise(act, sigma = 5)
#' \dontshow{rm(d,nscan,act,out)}
#' 
#' @keywords low-level noise

tasknoise <- function(act.image, sigma, type=c("gaussian","rician"), vee=1, verbose = T){
  
  if(missing(act.image)){
    stop("An activation array is required.")
  }
  if(missing(type)){
    type <- "gaussian"
  }
  
  if (verbose) cat("Task noise type: ", type, "\n")
  
  dim <- dim(act.image)
  if(length(dim)>4){
    stop("The activation array has more than 4 dimensions")
  }
  
  
  if(length(dim)==0){
    if(type=="gaussian"){
      noise <-rnorm(length(act.image), 0, sigma)
    } else {
      noise <- rrice(length(act.image), vee, sigma)
    }
  } else {
    if(type=="gaussian"){
      noise <- array(rnorm(prod(dim), 0, sigma), dim=dim)
    } else {
      noise <- array(rrice(prod(dim), vee, sigma), dim=dim)
    }
  }
  ix <- which(zapsmall(act.image)!=0)
  noise[-ix] <- 0
  noise <- noise*sigma/sd(noise)
  
  return(noise)
}
