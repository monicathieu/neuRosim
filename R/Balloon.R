Balloon <-
function(stim, T, acc){

	# stim: vector representing the presence/absence of a stimulus/activation
	# T: total duration of stimulus vector in seconds
	# acc: microtime resolution of stimulus vector in seconds

	require(deSolve, quiet=TRUE)
	
	kappa <- 2
	tau1 <- 3
        l <- -0.08
        rho <- 0.5
        k <- 3
        tauf <- 4
        tauh <- 0.242*tauf
        f1 <- 1.5
        deltat <- 1
        deltatf <- 0
        deltatm <- deltat - deltatf
        n <- 3
        m1 <- (f1-1)/n +1
        E0 <- 0.4
        v0 <- 0.03
        a1 <- 3.4
        a2 <- 1.0
	tauMTT <- 3
	tau <- 20
	alpha <- 0.4

	t <- seq(acc, T, acc)
	it <- c(1:(T/acc))

	par.I <- c(tau1, kappa)
	inhib <- function(t,y,p){
		yd1 <- (1/p[1])*(p[2]*stim[t] - (p[2]+1)*y[1])
		list(c(yd1))
	}
	I <- ode(c(0), it, inhib, par.I, method=rkMethod("ode45"))[,2]
	N <- stim - I
	F <- 1 + convolve((f1-1)*gammaHRF(t-deltatf,FWHM=4,verbose=FALSE),rev(N), type="o")
	M <- 1 + convolve((m1-1)*gammaHRF(t-deltatm,FWHM=4,verbose=FALSE),rev(N), type="o")
	E <- E0*M/F
	par.balloon <- c(E0, tauMTT, tau, alpha)
	de.balloon <- function(t,y,p){
		dv <- (F[t]-y[1]^(1/p[4]))/(p[2]+p[3])
		dq <- 1/p[2]*(F[t]*E[t]/p[1] - y[2]/y[1]*(y[1]^(1/p[4]) + p[3]/(p[2] + p[3])*(F[t] - y[1]^(1/p[4]))))
		list(c(dv,dq))
	}
	balloon <- ode(c(1,1),it,de.balloon,par.balloon,method=rkMethod("ode45"))
	V <- balloon[,2]
	Q <- balloon[,3]
	bold <- v0*(a1*(1-Q)-a2*(1-V))

	return(bold)
}

