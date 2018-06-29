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

gammaHRF <- function(x, FWHM = NULL, verbose = TRUE) {
  
  if (is.null(FWHM)) {
    if (verbose) warning("Default FWHM parameter used")
    FWHM <- 4
  }
  th <- 0.242*FWHM
  
  return (1/(th*factorial(3)) * (x/th)^3 * exp(-x/th))
}

balloon <- function(stim, totaltime, acc, par = NULL, verbose = TRUE){
    
    #require(deSolve, quietly=TRUE)
    
    if(missing(par)){
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
    inhib <- function(t,y,p){
      yd1 <- (1/p[1])*(p[2]*stim[t] - (p[2]+1)*y[1])
      list(c(yd1))
    }
    I <- ode(c(0), it, inhib, par.I, method=rkMethod("ode45"))[,2]
    N <- stim - I
    F <- 1 + convolve((par$f1-1)*gammaHRF(t-par$deltatf,FWHM=par$tauf,verbose=FALSE),rev(N), type="o")
    M <- 1 + convolve((par$m1-1)*gammaHRF(t-par$deltatm,FWHM=par$taum,verbose=FALSE),rev(N), type="o")
    E <- par$E0*M/F
    par.balloon <- c(par$E0, par$tauMTT, par$tau, par$alpha)
    de.balloon <- function(t,y,p){
      dv <- (F[t]-y[1]^(1/p[4]))/(p[2]+p[3])
      dq <- 1/p[2]*(F[t]*E[t]/p[1] - y[2]/y[1]*(y[1]^(1/p[4]) + p[3]/(p[2] + p[3])*(F[t] - y[1]^(1/p[4]))))
      list(c(dv,dq))
    }
    balloon <- ode(c(1,1),it,de.balloon,par.balloon,method=rkMethod("ode45"))
    V <- balloon[,2]
    Q <- balloon[,3]
    bold <- par$V0*(par$a1*(1-Q)-par$a2*(1-V))
    
    return(bold)
  }

