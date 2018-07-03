
#' Double-gamma Haemodynamic reponse function
#' Specifies a double-gamma variate haemodynamic response function for the given
#' time vector and parameters.
#' 
#' @export
#' 
#' @param x Time vector in seconds.
#' @param param List of parameters of the haemodynamic response function.
#' The list should contain the following:
#' \describe{
#'   \item{a1}{Delay of response relative to onset (default: 6)}
#'   \item{a2}{Delay of undershoot relative to onset (default:12)}
#'   \item{b1}{Dispersion of response (default:0.9)}
#'   \item{b2}{Dispersion of undershoot (default:0.9)}
#'   \item{c}{Scale of undershoot (default:0.35)}
#'   }
#' @param verbose If \code{TRUE}, warnings are displayed.
#' @return Vector representing the values of the function for the given
#' time vector and parameters.
#' 
#' @references 
#' Friston, KJ, Fletcher, P, Josephs, O, Holmes, AP, Rugg, MD and Turner, R (1998). Event-related fMRI: Characterising differential responses. NeuroImage, 7, 30-40.
#' Glover, GH (1999). Deconvolution of impulse response in event-related BOLD fMRI. NeuroImage, 9, 416-429.
#' 
#' @seealso \code{\link{gammaHRF}}, \code{\link{balloon}}
#' 
#' @examples
#' t <- 1:100
#' out <- canonicalHRF(t, verbose=FALSE)
#' \dontshow{rm(out, t)}
#' 
#' @keywords low-level activation

canonicalHRF <- function(x, param = NULL, verbose = TRUE){
  if(is.null(param)){
    if(verbose) warning("Default parameters for HRF are used")
    param <- list(a1 = 6,
                  a2 = 12,
                  b1 = 0.9,
                  b2 = 0.9,
                  c = 0.35)
  }
  d1 <- param$a1*param$b1		# time to peak of response
  d2 <- param$a2*param$b2		# time to peak of undershoot
  
  return ((x/d1)^param$a1*exp(-(x-d1)/param$b1) - param$c*(x/d2)^param$a2*exp(-(x-d2)/param$b2))
}

#' Single Gamma Haemodynamic response function
#' 
#' Specifies a Gamma variate haemodynamic response function for the given time vector and FWHM.
#' 
#' @export
#' 
#' @param x Time vector in seconds.
#' @param FWHM Full Width Half Maximum of the Gamma variate function.
#' @param verbose If \code{TRUE}, warnings are displayed.
#' @return Vector representing the values of the function for the given time vector and FWHM.
#' 
#' @references Buxton, RB, Uludag, K, Dubowitz, DJ and Liu, TT (2004). Modeling the hemodynamic response to brain activation. NeuroImage, 23, S220-S233.
#' 
#' @seealso \code{\link{canonicalHRF}}, \code{\link{balloon}}
#' 
#' @examples
#' t <- 1:100
#' out <- gammaHRF(t, verbose=FALSE)
#' \dontshow{rm(t,out)}
#' 
#' @keywords low-level activation

gammaHRF <- function(x, FWHM = NULL, verbose = TRUE) {
  
  if (is.null(FWHM)) {
    if (verbose) warning("Default FWHM parameter used")
    FWHM <- 4
  }
  th <- 0.242*FWHM
  
  return (1/(th*factorial(3)) * (x/th)^3 * exp(-x/th))
}

#' Balloon model
#' 
#' Generates the BOLD signal based on the Balloon model of Buxton et al. (2004).
#' 
#' @export
#' @import deSolve
#' @importFrom stats convolve
#' 
#' @param stim Vector representing the presence/absence (1-0 coding) of a stimulus/activation in seconds.
#' @param totaltime Total duration of stimulus vector in seconds.
#' @param acc Microtime resolution of stimulus vector in seconds.
#' @param par List representing the parameters of the Balloon model. The list should contain the following:
#' \describe{
#'   \item{kappa}{Inhibitory gain factor}
#'   \item{tau1}{Inhibitory time constant}
#'   \item{tauf}{FWHM of CBF impulse response}
#'   \item{taum}{FWHM of CMRO2 impulse resonse}
#'   \item{deltat}{Delay of CBF relative to CMRO2 response}
#'   \item{n}{Steady-state flow metabolism relation}
#'   \item{f1}{Normalized CBF response to sustained neural activation}
#'   \item{tauMTT}{Transit time through the balloon}
#'   \item{tau}{Viscoelastic time constant}
#'   \item{alpha}{Steady-state flow-volume relation}
#'   \item{E0}{baseline O2 extraction fraction}
#'   \item{V0}{baseline blood volume}
#'   \item{a1}{weight for deoxyHb change}
#'   \item{a2}{weight for blood volume change}
#'   }
#' @param verbose If \code{TRUE}, warnings are displayed.
#' @return Vector representing the values of the BOLD signal for the given stimulus
#' vector and microtime resolution.
#' 
#' @details Based on the provided stimulus boxcar function, a neural activation function
#' is generated that enters the Balloon model to generate a BOLD response. 
#' The microtime resolution ensures a high-precision generation of the response. More details
#' can be found in Buxton et al. (2004).
#' 
#' @references Buxton, RB, Uludag, K, Dubowitz, DJ and Liu, TT (2004). Modeling the hemodynamic response to brain activation. NeuroImage, 23, S220-S233.
#' 
#' @seealso \code{\link{canonicalHRF}}, \code{\link{gammaHRF}}
#' 
#' @examples
#' \dontrun{
#' s <- rep(rep(0,10), rep(1,10), 5)
#' T <- 100
#' it <- 0.1
#' out <- balloon(s, T, it) 
#' }
#' 
#' @keywords low-level activation

balloon <- function(stim, totaltime, acc, par = list(), verbose = TRUE){
    
    #require(deSolve, quietly=TRUE)
    
    if(length(par) == 0){
      if(verbose) warning("Default parameter values are used. See ?balloon.fnc for more details.")
    
      par <- list(kappa = 2,
                  tau1 = 3,
                  tauf = 4,
                  taum = 4,
                  f1 = 1.5,
                  deltat = 1,
                  n = 3,
                  E0 = 0.4,
                  V0 = 0.03,
                  a1 = 3.4,
                  a2 = 1.0,
                  tauMTT = 3,
                  tau = 20,
                  alpha = 0.4)
    }
    
    par$deltatf <- 0
    par$deltatm <- par$deltat - par$deltatf
    par$m1 <- (par$f1-1)/par$n +1
    
    t <- seq(acc, totaltime, acc)
    it <- c(1:(totaltime/acc))
    
    par.I <- c(par$tau1, par$kappa)

    I <- ode(c(0), it, inhib, par.I, method=rkMethod("ode45"))[,2]
    N <- stim - I
    F <- 1 + convolve((par$f1-1)*gammaHRF(t-par$deltatf,FWHM=par$tauf,verbose=FALSE),rev(N), type="o")
    M <- 1 + convolve((par$m1-1)*gammaHRF(t-par$deltatm,FWHM=par$taum,verbose=FALSE),rev(N), type="o")
    E <- par$E0*M/F
    par.balloon <- c(par$E0, par$tauMTT, par$tau, par$alpha)

    balloon <- ode(c(1,1),it,de.balloon,par.balloon,method=rkMethod("ode45"))
    V <- balloon[,2]
    Q <- balloon[,3]
    bold <- par$V0*(par$a1*(1-Q)-par$a2*(1-V))
    
    return(bold)
}

inhib <- function(t,y,p,stim){
  yd1 <- (1/p[1])*(p[2]*stim[t] - (p[2]+1)*y[1])
  list(c(yd1))
}

de.balloon <- function(t,y,p,E){
  dv <- (F[t]-y[1]^(1/p[4]))/(p[2]+p[3])
  dq <- 1/p[2]*(F[t]*E[t]/p[1] - y[2]/y[1]*(y[1]^(1/p[4]) + p[3]/(p[2] + p[3])*(F[t] - y[1]^(1/p[4]))))
  list(c(dv,dq))
}